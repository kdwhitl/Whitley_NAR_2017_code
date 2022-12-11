% Author: Kevin Whitley
% Date created: 111003
% Date last modified: 120709

% This function uses a T test to locate steps in a time trace.
% 
% Input:
% time - Array of all times in time trace
% allsteps - Array of all Y values in time trace
% pcutoff - Probability cutoff. Steps will only be counted if the
%   probability of their distribution is above this value. A probability
%   trace can be plotted to view this cutoff value.
% window - Window to use for T test to locate steps (in seconds)
%
% Output:
% steppedDown - Contains indexes of steps down
% steppedUp - Contains indexes of steps up
% p - p values

function [steppedDown, steppedUp, p] = ttesttrace2(time, allsteps, pcutoff, window)

% This part was / is for testing purposes.
if nargin == 0
    allsteps = [ones(500,1); zeros(100,1); ones(200,1); 2.*ones(300,1); 5.*ones(50,1); zeros(50,1); ones(200,1); zeros(200,1); ones(300,1)];
    allsteps = allsteps';
    allsteps = awgn(allsteps,10);
    
    time = 0:length(allsteps)-1;
    
    bw = 5;
    pcutoff = 3;
    window = 3; % Window size (sec)
else
    bw = round(1/(time(2)-time(1))); % Bandwidth
end

%% Find steps

windowind = floor(window*bw); % Window size (indexes)

allsteps = allsteps-mean(allsteps);
% allsteps = smooth(time,allsteps,smoother);
% allsteps = allsteps';

% Perform T test, get probabilities
p = ones(length(windowind+1:length(time)-windowind),1); % Initialize p
for iii = windowind+1:length(time)-windowind
    [~, p(iii)] = ttest2(allsteps(iii-windowind:iii),allsteps(iii:iii+windowind));
end

% Find middle of spike in probability. This corresponds to the approximate
% step locations.
meanind = 1;
approxstepLocs = zeros(1,length(time)-windowind); % Initialize
tt = 1;
while tt <= length(time)-windowind
    indAboveCut = 1;
    allAboveCut = [];
    while length(p) >= tt && -log10(p(tt)) > pcutoff
        allAboveCut(indAboveCut) = tt;
        tt = tt + 1;
        indAboveCut = indAboveCut + 1;
    end
    if ~isempty(allAboveCut)
        lowestP = min(p(allAboveCut));
        for dd = min(allAboveCut):max(allAboveCut)
            if p(dd)==lowestP
                approxstepLocs(meanind) = dd;
            end
        end
        meanind = meanind + 1;
    else
        tt = tt + 1;
    end
end
approxstepLocs(approxstepLocs==0) = []; % Remove zeros

% Identify whether the steps go up or down
downstep = 1;
upstep = 1;
steppedDown = zeros(length(approxstepLocs),1); % Initialize
steppedUp = zeros(length(approxstepLocs),1);
wiggleroom = ceil(bw/2); % Leave some points before and after step transition out of averaging
for bb = 1:length(approxstepLocs)
    if approxstepLocs(bb)-windowind <= 0
        beforestepind = 1;
    else
        beforestepind = approxstepLocs(bb)-windowind:approxstepLocs(bb)-wiggleroom;
    end
    if approxstepLocs(bb)+windowind > length(allsteps)
        afterstepind = length(allsteps);
    else
        afterstepind = approxstepLocs(bb)+wiggleroom:approxstepLocs(bb)+windowind;
    end
    
    stepsize1 = mean(allsteps(afterstepind)) - mean(allsteps(beforestepind));
    if stepsize1 <= 0
        steppedDown(downstep) = approxstepLocs(bb);
        downstep = downstep + 1;
    else
        steppedUp(upstep) = approxstepLocs(bb);
        upstep = upstep + 1;
    end
end
steppedDown(steppedDown==0) = []; % Remove zeros
steppedUp(steppedUp==0) = [];

%% Analyze steps down
% 
% stepDown = [];
% ptsbeforedown = [];
% ptsafterdown = [];
% for nnn = 1:length(steppedDown)
%     tdownevent = steppedDown(nnn);
%     
%     eventbefore = time(steppedDown(steppedDown < tdownevent));
%     if isempty(eventbefore)
%         eventbefore = time(1);
%     end
%     closestdownbefore = min(time(tdownevent) - eventbefore);
%     
%     eventbefore = time(steppedUp(steppedUp < tdownevent));
%     if isempty(eventbefore)
%         eventbefore = time(1);
%     end
%     closestupbefore = min(time(tdownevent) - eventbefore);
%     
%     mintbeforedown = min([closestdownbefore; closestupbefore]);
%     
%     availptsbefore = mintbeforedown * bw - 2;
%     useptsbefore = floor(min([availptsbefore; windowind]));
%     
%     trapokbefore = availptsbefore > bw;
%     
%     eventafter = time(steppedDown(steppedDown > tdownevent));
%     if isempty(eventafter)
%         eventafter = time(end);
%     end
%     closestdownafter = min(eventafter - time(tdownevent));
%     
%     eventafter = time(steppedUp(steppedUp > tdownevent));
%     if isempty(eventafter)
%         eventafter = time(end);
%     end
%     closestupafter = min(eventafter - time(tdownevent));
%     
%     mintafterdown = min([closestdownafter; closestupafter]);
%     
%     availptsafter = mintafterdown * bw - 2;
%     useptsafter = floor(min([availptsafter; windowind]));
%     
%     trapokafter = availptsafter > bw;
%     
%     if ~(trapokbefore && trapokafter)
%         continue % If there aren't enough points to average, skip this event
%     end
%     
%     t_before = time(time < time(tdownevent));
%     indexbefore = length(t_before)-useptsbefore:length(t_before)-1;
%     ptsbeforedown = [ptsbeforedown; t_before(indexbefore)' allsteps(indexbefore)'];
%     
%     t_after = time(time > time(tdownevent));
%     step_after = allsteps(time > time(tdownevent));
%     if useptsafter > 0
%         indexafter = (1:useptsafter)+1;
%     else
%         indexafter = 1;
%     end
%     if indexafter(end) > length(t_after) % Handles case at end of data
%         thend = indexafter(end) - length(t_after);
%         indexafter = indexafter(1:end-thend);
%     end
%     ptsafterdown = [ptsafterdown; t_after(indexafter)' step_after(indexafter)'];
%     
%     dextbefore = allsteps(1, time <= time(tdownevent));
%     dextafter = allsteps(1, time > time(tdownevent));
%     
%     dextbefore1 = dextbefore(indexbefore);
%     extmeanbefore = mean(dextbefore1);
%     dextafter1 = dextafter(indexafter);
%     extmeanafter = mean(dextafter1);
%     
%     deltaext = extmeanafter - extmeanbefore;
%     if ~isempty(steppedDown)
%         stepdownind = find(time==time(tdownevent),1,'last');
%     end
%     stepDown = [stepDown; stepdownind time(tdownevent) deltaext];
% end
% 
% %% Analyze steps up
% 
% stepUp = [];
% ptsbeforeup = [];
% ptsafterup = [];
% for pp = 1:length(steppedUp)
%     tupevent = steppedUp(pp);
%     
%     eventbefore = time(steppedUp(steppedUp < tupevent));
%     if isempty(eventbefore)
%         eventbefore = time(1);
%     end
%     closestupbefore = min(time(tupevent) - eventbefore);
%     
%     eventbefore = time(steppedDown(steppedDown < tupevent));
%     if isempty(eventbefore)
%         eventbefore = time(1);
%     end
%     closestdownbefore = min(time(tupevent) - eventbefore);
%     
%     mintbeforeup = min([closestupbefore; closestdownbefore]);
%     
%     availptsbefore = mintbeforeup * bw - 2;
%     useptsbefore = floor(min([availptsbefore; windowind]));
%     
%     trapokbefore = availptsbefore > bw;
%     
%     eventafter = time(steppedUp(steppedUp > tupevent));
%     if isempty(eventafter)
%         eventafter = time(end);
%     end
%     closestupafter = min(eventafter - time(tupevent));
%     
%     eventafter = time(steppedDown(steppedDown > tupevent));
%     if isempty(eventafter)
%         eventafter = time(end);
%     end
%     closestdownafter = min(eventafter - time(tupevent));
%     
%     mintafterup = min([closestupafter; closestdownafter]);
%     
%     availptsafter = mintafterup * bw - 2;
%     useptsafter = floor(min([availptsafter; windowind]));
%     
%     trapokafter = availptsafter > bw;
%     
%     if ~(trapokbefore && trapokafter)
%         continue % If there aren't enough points to average, skip this event
%     end
%     
%     t_before = time(time < time(tupevent));
%     indexbefore = length(t_before)-useptsbefore:length(t_before)-1;
%     ptsbeforeup = [ptsbeforeup; t_before(indexbefore)' allsteps(indexbefore)'];
%     
%     t_after = time(time > time(tupevent));
%     step_after = allsteps(time > time(tupevent));
%     if useptsafter > 0
%         indexafter = (1:useptsafter)+1;
%     else
%         indexafter = 1;
%     end
%     if indexafter(end) > length(t_after) % Handles case at end of data
%         thend = indexafter(end) - length(t_after);
%         indexafter = indexafter(1:end-thend);
%     end
%     ptsafterup = [ptsafterup; t_after(indexafter)' step_after(indexafter)'];
%     
%     dextbefore = allsteps(1, time <= time(tupevent));
%     dextafter = allsteps(1, time > time(tupevent));
%     
%     dextbefore1 = dextbefore(indexbefore);
%     extmeanbefore = mean(dextbefore1);
%     dextafter1 = dextafter(indexafter);
%     extmeanafter = mean(dextafter1);
%     
%     deltaext1 = extmeanafter - extmeanbefore;
%     if ~isempty(steppedUp)
%         stepupind = find(time==time(tupevent),1,'last');
%     end
%     stepUp = [stepUp; stepupind time(tupevent) deltaext1];
% end

% Eliminate spikes in probability trace
% newUp = [];
% newDown = [];
% newupind = 1;
% newdownind = 1;
% for mmm = 1:length(stepDown)
%     nextUp = min(stepUp(stepUp >= stepDown(mmm)));
%     nextDown = min(stepDown(stepDown > stepDown(mmm)));
%     if nextDown < nextUp
%         continue
%     end
%     if ~isempty(nextUp) && abs(stepDown(mmm)-nextUp)>lifetimecutoff
%         if newupind~=1 && nextUp~=newUp(newupind-1)
%             newUp(newupind) = nextUp;
%             newupind = newupind + 1;
%         elseif ~isempty(nextUp) && newupind == 1
%             newUp(newupind) = nextUp;
%             newupind = newupind + 1;
%         end
%         if ~isempty(nextDown) && newdownind~=1 && nextDown~=newDown(newdownind-1)
%             newDown(newdownind) = nextDown;
%             newdownind = newdownind + 1;
%         elseif ~isempty(nextDown) && newdownind == 1
%             newDown(newdownind) = nextDown;
%             newdownind = newdownind + 1;
%         end
%     end
% end

%% Plot results

if 0
    figSize = [25 49 900 600];
    figure('NumberTitle', 'off', 'Name', 'Ttest', 'Pos', figSize)

    subplot(2,1,1)
    plot(time,allsteps,'k')
    
    if ~isempty(steppedDown)
        hold on
        down = plot(time(steppedDown), max(allsteps)+max(allsteps)*0.5, 'xg', 'MarkerSize', 10, 'linew', 2);
%         hold on
%         plot(ptsbeforedown(:,1),ptsbeforedown(:,2),'^m','MarkerSize',3)
%         hold on
%         plot(ptsafterdown(:,1),ptsafterdown(:,2),'or','MarkerSize',3)
    end
    if ~isempty(steppedUp)
        hold on
        up = plot(time(steppedUp), max(allsteps)+max(allsteps)*0.5, '+g', 'MarkerSize', 10, 'linew', 2);
%         hold on
%         plot(ptsbeforeup(:,1),ptsbeforeup(:,2),'^g','MarkerSize',3)
%         hold on
%         plot(ptsafterup(:,1),ptsafterup(:,2),'ob','MarkerSize',3)
    end
    xlim([time(1) time(end)])
    ax1 = gca;
    ylim([min(allsteps)+min(allsteps)*0.5 max(allsteps-mean(allsteps))+max(allsteps-mean(allsteps))*1.5])
    ylabel('\DeltaExtension')
    
    subplot(2,1,2)
    plot(time(1:end-windowind),-log10(p), 'k')
    
    ax2 = gca;
    ylim([min(-log10(p)) max(-log10(p))])
    hold on
    cutoffmatrix = zeros(length(time),1);
    cutoffmatrix(cutoffmatrix==0) = pcutoff;
    
    plot(time, cutoffmatrix, '--r')
    ylabel('Probability of non-zero mean')
    
%     legend([up(1), down(1)], 'Steps up', 'Steps down')
    
    linkaxes([ax1 ax2], 'x')
end