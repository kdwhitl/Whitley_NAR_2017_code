function [avediff, upsteps, downsteps] = findstepsK(data,avewindow,smthx,transthresh)
%find up/down steps in data
%averages and produces running derivative
%then finds threholded local min/max of a 1D array using matlab's 2D image min/max routines
%100420 mjc

l = length(data);

shmatahead = zeros(avewindow,l);
shmatbehind = zeros(avewindow,l);

% Make several columns of data, each one shifted by 1 from the last.
for awi = 1:avewindow
    shmatahead(awi,:) = circshift(data',-awi)';
    shmatbehind(awi,:) = circshift(data',awi-1)';
end

% Take derivative. This assumes a constant spacing between points.
aheadmean = sum(shmatahead)/avewindow; % avewindow is [indexes].
behindmean = sum(shmatbehind)/avewindow;

avediff = aheadmean - behindmean;

%smooth the derivative a bit if necessary
if smthx ~= 0
    avediff = smooth(avediff,smthx,'moving')';
end

%now find local min/max

% Find middle of spike in probability. This corresponds to the approximate
% step locations. Not very elegant, but it works. This seems very ugly with
% all the while loops, but it finds events better than the other method of
% searching the whole trace for local max/min (what we want is the GLOBAL
% maximum for EACH threshold crossing). KW 141104
meanind = 1;
tt = 1;
steppedUp = [];
while tt <= length(data)-avewindow
    indAboveCut = 1;
    allAboveCut = [];
    while length(avediff) >= tt && avediff(tt) > transthresh % Steps up
        allAboveCut(indAboveCut) = tt;
        tt = tt + 1;
        indAboveCut = indAboveCut + 1;
    end
    if ~isempty(allAboveCut)
        avedifflocal = zeros(1,length(avediff));
        avedifflocal(allAboveCut) = avediff(allAboveCut);
        [~, steppedUp(meanind)] = max(avedifflocal);
        meanind = meanind + 1;
    else
        tt = tt + 1;
    end
end
steppedUp(steppedUp==0) = [];

tt2 = 1;
meanind2 = 1;
steppedDown = [];
while tt2 <= length(data)-avewindow
    indBelowCut = 1;
    allBelowCut = [];
    while length(avediff) >= tt2 && avediff(tt2) < -transthresh % Steps down
        allBelowCut(indBelowCut) = tt2;
        tt2 = tt2 + 1;
        indBelowCut = indBelowCut + 1;
    end
    if ~isempty(allBelowCut)
        avedifflocal = zeros(1,length(avediff));
        avedifflocal(allBelowCut) = avediff(allBelowCut);
        [~, steppedDown(meanind2)] = min(avedifflocal);
        meanind2 = meanind2 + 1;
    else
        tt2 = tt2 + 1;
    end
end
steppedDown(steppedDown==0) = [];

imupsteps = zeros(1,length(data));
imupsteps(steppedUp) = 1;

imdownsteps = zeros(1,length(data));
imdownsteps(steppedDown) = 1;

% im = repmat(avediff,1,1);
% imupstepscut = (im > transthresh).*im;
% imdownstepscut = (im < -transthresh).*im;
% 
% imupsteps = imregionalmax(imupstepscut);
% imdownsteps = imregionalmin(imdownstepscut);

% if max(imupstepscut(:)) == 0
%     upsteps = imupsteps(1,:)*0;
% else
%     upsteps = imupsteps(1,:);
% end
% 
% if min(imdownstepscut(:)) == 0
%     downsteps = imdownsteps(1,:)*0;
% else
%     downsteps = imdownsteps(1,:);
% end

if max(imupsteps(:)) == 0
    upsteps = imupsteps(1,:)*0;
else
    upsteps = imupsteps(1,:);
end

if min(-imdownsteps(:)) == 0
    downsteps = imdownsteps(1,:)*0;
else
    downsteps = imdownsteps(1,:);
end

%filter false beginning and ending events
mask = ones(1,length(imupsteps));
mask(1:avewindow) = zeros(1,avewindow);
mask(end-(avewindow-1):end) = zeros(1,avewindow);

upsteps = upsteps .* mask;
downsteps = downsteps .* mask;
