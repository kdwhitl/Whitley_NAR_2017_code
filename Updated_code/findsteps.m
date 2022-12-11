function [avediff, upsteps, downsteps] = findsteps(data,avewindow,smthx,transthresh)
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

aheadmean = sum(shmatahead)/avewindow;
behindmean = sum(shmatbehind)/avewindow;

avediff = aheadmean - behindmean;

%smooth the derviative a bit if necessary
if smthx ~= 0
    avediff = smooth(avediff,smthx,'moving')';
end

%now find local min/max
im = repmat(avediff,1,1);
imupstepscut = (im > transthresh).*im;
imdownstepscut = (im < -transthresh).*im;

imupsteps = imregionalmax(imupstepscut);
imdownsteps = imregionalmin(imdownstepscut);

if max(imupstepscut(:)) == 0
    upsteps = imupsteps(1,:)*0;
else
    upsteps = imupsteps(1,:);
end

if min(imdownstepscut(:)) == 0
    downsteps = imdownsteps(1,:)*0;
else
    downsteps = imdownsteps(1,:);
end

%filter false beginning and ending events
mask = ones(1,length(upsteps));
mask(1:avewindow) = zeros(1,avewindow);
mask(end-(avewindow-1):end) = zeros(1,avewindow);

upsteps = upsteps .* mask;
downsteps = downsteps .* mask;
