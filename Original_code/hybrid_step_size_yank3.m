% Author: Matt Comstock
% Date created: 100114 (yymmdd)
% Date last modified: 150317 by Kevin Whitley

% For step-finding data: order in DataSets is 'Date', Data file, Cal file,
% [Time interval], Fluor threshold (Hz), Smoothing factor, APD averaging
% factor, T test threshold, [unlabeled binding and unbinding events], oligo
% concentration.
tic
global datadirectories
global beaddiameters
SaveFigDir = 'C:\Users\Kevin\Documents\UHR_Fleezer_Data\Oligo_data\';

%% Data sets
probe = 'seq6'; % NOTE: FLUORESCENCE THRESHOLD CHANGED 150317. THEN CHANGED BACK 150325.
assignin('base', 'probe', probe)
DataSets = hybridDataFile3(probe);

%% Initialize

% DSetsToAnalyze = 5;
DSetsToAnalyze = 1:length(DataSets);

overview = 0;
fluor_plots = 1; % Plot fluorescence
removeyankoffsets = 1;

% Use t test for step finding
step_find = 1; % Use T test to find events from trap trace. If OT = 1, this defaults to 1. Also, disables finding steps via fluorescence.
step_plots = 0; % Plot step-finding?

% Initialize structure variables
AllResults.binding = [];
AllResults.unbinding = [];
AllResults.lifetimes = [];
AllResults.rupture = [];
YankResults = [];

desired_trap_bw = 5; % [Hz], APD will be averaged to match
if strcmp(probe,'7merRepeatHairpin')
    desired_trap_bw = 10;
end

% avewindow = desired_trap_bw*3; % Fl averaging window, number of data points
% t_avewindow = 0.15; % [s] Fl averaging window.
% trapdpoints = desired_trap_bw*1.5; % Trap data points to include before/after hybrid event, 9 is 3 s
XAODCal = 123; % [nm/MHz]
avgptcutoff = desired_trap_bw * 0.5; % This basically determines how short an event can be for taking an extension change.
avplotfact = 20; % 'Downsampling' factor for plotting

% For old trap data
OT = 0; % Not using Old Trap unless otherwise specified. Only used for some 12mer data.
fixedMirrorY = 3.4; %[V]
convY = 505.5; % [nm/V]

Cy3 = [0 0.7 0]; % Green color

for ii = DSetsToAnalyze
%% Read data files

    Date =                  DataSets{ii}{1};
    datafnumber =           DataSets{ii}{2};
    calfilenumber =         DataSets{ii}{3};
    plotinterval =          DataSets{ii}{4};
    
    % This is now /100 because the APD signal is no longer multiplied by
    % BW. Please be careful of this! 150317. Changed back 150325.
    fltransthresh =         DataSets{ii}{5}; % [Hz]
    smthx =                 DataSets{ii}{6}; % Fl deriv smoothing factor
    
    if length(DataSets{ii})>=7
        APDAveFactor =      DataSets{ii}{7};
    else
        APDAveFactor = 10;
    end
    
    % Probability cutoff for t test (-log10(p value))
    if length(DataSets{ii})>=8
        probcutoff =        DataSets{ii}{8};
    else
        probcutoff = 4;
    end
    
    % List of unlabeled binding/unbinding events for T test
    if length(DataSets{ii})>=9
        unlabeled =         DataSets{ii}{9};
    else
        unlabeled = [];
    end
    
    % If [oligo] is not 10 nM
    if length(DataSets{ii}) >= 10
        oligo_conc =        DataSets{ii}{10};
    else
        oligo_conc = 10; % [nM] Concentration of probe for calculating kon.
    end
    
    % Window size for t test
    if length(DataSets{ii}) >= 11
        windowSize =        DataSets{ii}{11};
    else
        windowSize = 3; % [s]
    end
    
    % NE filter setting
    if length(DataSets{ii}) >= 12
        NE =                DataSets{ii}{12};
    else
        NE = 20;
    end
    
    % Green FB QPD setpoint
    if length(DataSets{ii}) >= 13
        FBSP =              DataSets{ii}{13};
    else
        FBSP = 90; % [mV]
    end
    
    % A few exceptions...
    if strcmp(probe,'12mer1') && ii>=110 && ii<=127 % Part of the 12mer1 data was saved with the old trap.
        OT = 1;
    elseif strcmp(probe,'12mer1') && (ii<110 || ii>127)
        OT = 0;
    end
    
    if strcmp(probe,'9merIntCy3\_NoSpacer') && datafnumber>=259 || strcmp(probe,'NoSpacer0mg')
%         step_find = 1;
    end
    
    % Get data
    if ~OT % Fleezers
        calrootfile = [Date '_' num2str(calfilenumber,'%03d') '.dat'];
        startpath = [datadirectories Date '\'];
        rootfile = [Date '_' num2str(datafnumber, '%03d') '.dat'];
        data = ReadMattFile_Wrapper(startpath, rootfile);
    else % Old trap
        calrootfile = [Date 'N' num2str(calfilenumber, '%02d') '.dat'];
        samplerate = 100;
        startpath = [rawDataPath Date '\'];
        rootfile = [Date 'N' num2str(datafnumber, '%02d') '.dat'];
        data = ReadFileK(samplerate, startpath, rootfile);
        step_find = 1;  % If using old trap data, you must be using the t test, not fluorescence.
    end
    calfilename = ['cal' calrootfile(1:(end-4))];
    
    % If the calibration is already in memory, don't have to redo
    if str2double(Date) < 140821
        beaddiameters = [790 880];
    else
        beaddiameters = [950 880];
    end
    if ~exist(calfilename, 'var') && ~OT % Then do calibration
        CalibrateStrippedWrapper(Date, calrootfile, 0, 0);
    elseif ~exist(calfilename, 'var') && OT
        CalibrateStrippedWrapperOT(calrootfile, 0, 0, Date);
    end
    cal = evalin('base', calfilename);
    evalin('base','clear([''fit'' calrootfile(1:(end-4))])') % Clear out this (unused) structure variable to save memory.
    
    if OT % Account for different pulling directions with old trap and fleezers data.
        cal.kappaA = cal.kappaAY;
        cal.kappaB = cal.kappaBY;
        cal.alphaA = cal.alphaAY;
        cal.alphaB = cal.alphaBY;
        cal.Aoffset = cal.AYoffset;
        cal.Boffset = cal.BYoffset;
        data.A = data.A_Y;
        data.B = data.B_Y;
    else
        cal.kappaA = cal.kappaAX;
        cal.kappaB = cal.kappaBX;
        cal.alphaA = cal.alphaAX;
        cal.alphaB = cal.alphaBX;
        cal.Aoffset = cal.AXoffset;
        cal.Boffset = cal.BXoffset;
        data.A = data.A_X;
        data.B = data.B_X;
    end

%% Integrate/downsample fluorescence data

    if ~OT % Fleezers data
        if str2double(Date) < 160101 % Changed to APD2 after realignment
            apd = apd_integrate(data.apd1, APDAveFactor);
        else
            apd = apd_integrate(data.apd2, APDAveFactor);
        end
        t_apd = downsample(data.apdtime, APDAveFactor);
        
        t_apdav = downsample(data.apdtime, avplotfact); % For plotting
        if str2double(Date) < 160101 % Changed to APD2 after realignment
            apdav = apd_integrate(data.apd1, avplotfact); % [Hz]
        else
            apdav = apd_integrate(data.apd2, avplotfact);
        end
        apdavBW = 1/(t_apdav(2)-t_apdav(1)); % [Hz]
        apdav = apdav * apdavBW / 1000; % [kHz]
        
        ApdBW = 1/(t_apd(2)-t_apd(1)); % [Hz]
        apdr = apd * ApdBW; % [Hz]
        apdr1 = apdr./1000; % [kHz]
    end
        
%% Offset removal part 1: Determine what the offset correction is at each force

    % Small exception. No offset files in the first 94 time traces of this
    % data set.
    if strcmp(probe,'9merCy3')
        if ii >= 95
            removeyankoffsets = 1;
        else
            removeyankoffsets = 0;
        end
    end
    
    if removeyankoffsets
        offsetfilenumber = calfilenumber + 1;
        if ~OT
            offsetrootfile = [Date '_' num2str(offsetfilenumber,'%03d') '.dat'];
        else
            offsetrootfile = [Date 'N' num2str(offsetfilenumber,'%02d') '.dat'];
        end
        offsetfilename = ['offset' offsetrootfile(1:(end-4))];
        
        % If the offset is already in memory, don't redo
        if ~exist(offsetfilename, 'var') % Then compute offset
            if ~OT
                offsetdata = ReadMattFile_Wrapper(startpath, offsetrootfile);
            else % Old trap
                offsetdata = ReadFileK(samplerate, startpath, rootfile);
            end
            if isfield(offsetdata,'nscanfb')
                offsetdata = separate_raster_data(offsetdata);
            end
            
            if ~OT
                % Take mean of multiple scans (uncalibrated)
                offset.axoffraw = mean(offsetdata.A_X);
                offset.ayoffraw = mean(offsetdata.A_Y);
                offset.bxoffraw = mean(offsetdata.B_X);
                offset.byoffraw = mean(offsetdata.B_Y);
                offset.nu = (offsetdata.t2x-offsetdata.t1x) - offsetdata.nu;
            else
                offset.axoffraw = offsetdata.A_X;
                offset.ayoffraw = offsetdata.A_Y;
                offset.bxoffraw = offsetdata.B_X;
                offset.byoffraw = offsetdata.B_Y;
                % offset.nu = (offsetdata.Mirror_Y);
                offset.nu = [];
            end
            assignin('base', offsetfilename, offset);
        else
            offset = evalin('base', offsetfilename);
        end
        
        % Correlate trap separations [MHz] to tensions [pN] (rough estimate)
        if ~OT
            yanknu = [14.6 14.95 15.35 15.6 15.8 16.15 16.4 16.65]; % [MHz]
            yankT = [3 5 8 10 12 15 17.5 20]; % [pN]
        else
            yanknu = [6.945 6.98 7.15]; % [MHz]
            yankT = [15 20 25]; % [pN]
        end
        
        % Get set of offsets for each trap separation spot (yanknu)
        yankoffax = [];
        yankoffbx = [];
        yankoffay = [];
        yankoffby = [];
        for ynu = yanknu
            [~, ind] = min(abs(ynu - offset.nu)); % Find trap separations in offset file corresponding to those listed above
            if ~OT
                yankoffax = [yankoffax offset.axoffraw(ind)];
                yankoffbx = [yankoffbx offset.bxoffraw(ind)];
            else
                yankoffay = [yankoffay offset.ayoffraw(ind)];
                yankoffby = [yankoffby offset.byoffraw(ind)];
            end
        end
    else
        yanknu = [];
    end
    
    % Get a calibration for each trap separation spot
    APos = zeros(length(yanknu)+1, length(data.A));
    BPos = APos;
    PDiff = APos;
    
    % No offset correction
    APos(1,:) = (data.A-cal.Aoffset) * cal.alphaA;
    BPos(1,:) = (data.B-cal.Boffset) * cal.alphaB;
    
    % Get difference in bead / trap positions
    if FB && ~OT % Force feedback
        PDiff(1,:) = -data.trappos1*XAODCal;
    elseif FB && OT
        PDiff(1,:) = (data.Mirror_Y-fixedMirrorY)*convY; % This should not be negative. Check DNAStretch2.m for reference. 120619
    elseif ~FB
        PDiff(1,:) = BPos(1,:) - APos(1,:);
    end
    
%% Offset removal part 2: Actually remove offsets from data

    TrapAveFactor = round((1 / data.sampperiod) / desired_trap_bw);
    t_trap = downsample(data.time, TrapAveFactor);
    traponly = downsample(data.trappos1, TrapAveFactor);
    bw_trap = 1 / (t_trap(2) - t_trap(1));
    
    Dext = downsample(PDiff(1,:), TrapAveFactor);
    [~, c] = size(Dext);
    dext = zeros(length(yanknu)+1, c);
    Aforce = dext;
    Bforce = dext;
    
    dext(1,:) = Dext;
    Aforce(1,:) = downsample(APos(1,:)*cal.kappaA, TrapAveFactor);
    Bforce(1,:) = downsample(BPos(1,:)*cal.kappaB, TrapAveFactor);
    
    % Yanks corrected for offsets
    if removeyankoffsets
        for ppp = 1:length(yanknu)
            APos(ppp+1,:) = (data.A-yankoffax(ppp)) * cal.alphaA;
            BPos(ppp+1,:) = (data.B-yankoffbx(ppp)) * cal.alphaB;
            Aforce(ppp+1,:) = downsample(APos(ppp+1,:)*cal.kappaA, TrapAveFactor);
            Bforce(ppp+1,:) = downsample(BPos(ppp+1,:)*cal.kappaB, TrapAveFactor);
            PDiff(ppp+1,:) = BPos(ppp+1,:) - APos(ppp+1,:);
            dext(ppp+1,:) = downsample(PDiff(ppp+1,:), TrapAveFactor); % Each row in dext is now change in extension with a different offset subtracted.
        end
    end
    
    clear APos BPos PDiff data % Save memory
    
    if ~OT
        avgF = (Aforce(1,:)-Bforce(1,:))/2;
    else
        avgF = (Bforce(1,:)-Aforce(1,:))/2;
    end
    
%% Downsample trap data

    % For step finding code (and later on) we just want a small segment of
    % data, so cut the data down to plotinterval.
    t_local = t_trap(t_trap>plotinterval(1) & t_trap<plotinterval(2));
    dextLocal = dext(1,t_trap>plotinterval(1) & t_trap<plotinterval(2));
    
    t_apdlocal = t_apd(t_apd>plotinterval(1) & t_apd<plotinterval(2));
    apdlocal = apdr(t_apd>plotinterval(1) & t_apd<plotinterval(2));
    
    AF_local = Aforce(:,t_trap>plotinterval(1) & t_trap<plotinterval(2));
    BF_local = Bforce(:,t_trap>plotinterval(1) & t_trap<plotinterval(2));
        
    F_local = (AF_local(1,:)-BF_local(1,:))/2;

%% Plot just simple fluorescence and trap data (no analysis)

    if (fluor_plots || step_plots) && overview % Fleezers data
        
        figure('NumberTitle', 'off', 'FileName', [SaveFigDir Date '_' num2str(datafnumber, '%03d')...
            '_' probe], 'Name', [Date '_' num2str(datafnumber, '%03d') '_' probe ' Overview'])
        
        if fluor_plots
            subplot(3,1,1,'Ylim',[0 max(apdav)+max(apdav)*0.5]) %#ok<*UNRCH>
            box on
            
            plot(t_apdav, apdav, 'Color', Cy3)
            ylabel('Photon Rate (kHz)')
            ax = gca;
            
            apdleg = legend([num2str(apdavBW, '%10.1f') ' Hz']);
            title([probe '  ' Date '\_' num2str(datafnumber)])
            set(gca,'XTickLabel',[])
        end
        
        if fluor_plots
            subplot(312)
            box on
        elseif step_plots
            subplot(211)
            title([probe '  ' Date '\_' num2str(datafnumber)])
        end

        plot(t_trap, dext(1,:)-mean(dext(1,:)),'k')
        ylim([min(dext(1,:))-mean(dext(1,:)) max(dext(1,:))-mean(dext(1,:))])
        ylabel('\DeltaExtension (nm)')
        ax = [ax gca];
        dextleg = legend([num2str(bw_trap, '%10.1f') ' Hz']);
        set(gca,'XTickLabel',[])
        
        if fluor_plots
            subplot(3,1,3,'YLim',[min([Aforce(1,:) -Bforce(1,:)]) max([Aforce(1,:) -Bforce(1,:)])])
            box on
        elseif step_plots
            subplot(2,1,2,'YLim',[min([-Aforce(1,:) Bforce(1,:)]) max([-Aforce(1,:) Bforce(1,:)])])
        end
        
        if ~OT
            plot(t_trap,Aforce(1,:),'b',t_trap,-Bforce(1,:),'r',t_trap,avgF,'k')
            box on
        else
            plot(t_trap,-Aforce(1,:),'b',t_trap,Bforce(1,:),'r',t_trap,avgF,'k')
        end
        ylabel('Force (pN)')
        xlabel('Time (s)')
        forceleg = legend('F_A','F_B','AvgF');
        ax = [ax gca];
        
        linkaxes(ax, 'x')
    end
        
%% Find binding/unbinding events
        
    % Find events by taking derivative of fluorescence
    if ~OT
        % Find binding/unbinding steps in fl data
        avewindow = floor(ApdBW * t_avewindow); % [index]
        [apdravediff, bindingevents, unbindingevents] = findstepsK(apdlocal, avewindow, smthx, fltransthresh);
        
%         fpcutoff = 4;
%         fwinsize = 1;
%         [fluorDown, fluorUp, fprob] = ttesttrace2(t_apdlocal, apdlocal, fpcutoff, fwinsize);

        t_apdevents = t_apdlocal + (1/ApdBW * 1/2); % Now event times are set in the middle of the transition
        
        % These are now the times when binding/unbinding ocurred
%         t_binding = fluorDown ./ ApdBW;
%         t_unbinding = fluorUp ./ ApdBW;
        t_binding = bindingevents .* t_apdevents;
        t_unbinding = unbindingevents .* t_apdevents;
        
        % Get rid of zeros, then have list of times
        t_binding = t_binding(t_binding ~= 0);
        t_unbinding = t_unbinding(t_unbinding ~= 0);
        
        % Only keep events that are within the plot window
        t_binding = t_binding(t_binding > plotinterval(1) & t_binding < plotinterval(2));
        t_unbinding = t_unbinding(t_unbinding > plotinterval(1) & t_unbinding < plotinterval(2));
        
        % Numbers of events
        nbinding = length(t_binding);
        nunbinding = length(t_unbinding);
    end

    % Find events by paired t test of trap signal
    if step_find
        
        [stepDown, stepUp, prob] = ttesttrace2(t_local, dextLocal, probcutoff, windowSize);
        windowind = floor(windowSize*desired_trap_bw); % Window size (indexes)
        
        t_binding = t_trap(stepDown)+plotinterval(1);
        t_unbinding = t_trap(stepUp)+plotinterval(1);
        t_all = sort([t_binding t_unbinding]);
        
        nbinding = length(stepDown);
        nunbinding = length(stepUp);
    end
    
%% Find yanks
    
    traptransthresh = 1; % [pN]
    trapsmthx = 9; % Trap deriv smoothing factor
    yank_avewindow = 15; % [index]
    [~, yankingevents, unyankingevents] = findsteps(Aforce(1,:), yank_avewindow, trapsmthx, traptransthresh);
    
    t_yankingevents = t_trap + (1/bw_trap * 1/2); % Now event times are set in the middle of the transition
    
    % These are now the times when yanking/unyanking ocurred
    t_yanking = yankingevents .* t_yankingevents;
    t_unyanking = unyankingevents .* t_yankingevents;
    
    % Get rid of zeros, then have list of times
    t_yanking = t_yanking(t_yanking ~= 0);
    t_unyanking = t_unyanking(t_unyanking ~= 0);
    allyanks = [0 sort([t_yanking t_unyanking]) t_trap(end)]; % Times where all force changes happened (for ENTIRE trace, not just in plot window).
    
    % Only keep events that are within the plot window
    t_yanking = t_yanking(t_yanking>plotinterval(1) & t_yanking<plotinterval(2));
    t_unyanking = t_unyanking(t_unyanking>plotinterval(1) & t_unyanking<plotinterval(2));
    
    t_yanking = [plotinterval(1)+1 t_yanking]; % Should have at least one yanking event in here for offset removal. KW 120629
    
    % Identify size of each yank
    ypriorT = [];
    ymeasT = [];
    yapproxT = [];
    ynu = [];
    yankid = [];
    if removeyankoffsets
        for yk = 1:length(t_yanking) % Should at least run once to remove offset.
            [~, ind] = min(abs(t_yanking(yk) - t_trap));
            %             ypriorT = [ypriorT AXforce(1,ind-2)];
            %             ymeasT = [ymeasT AXforce(1,ind+2)];
            ypriorT = [ypriorT avgF(1,ind-2)];
            ymeasT = [ymeasT avgF(1,ind+2)];
                        
            [~, ind2] = min(abs(ymeasT(yk) - yankT));
            yapproxT = [yapproxT yankT(ind2)];
            ynu = [ynu yanknu(ind2)];
            yankid = [yankid ind2];
        end
        for yk = 1:length(t_unyanking)
            [~, ind] = min(abs(t_unyanking(yk) - t_trap));
            %                 ypriorT = [ypriorT AXforce(1,ind-2)];
            %                 ymeasT = [ymeasT AXforce(1,ind+2)];
            ypriorT = [ypriorT avgF(1,ind-2)];
            ymeasT = [ymeasT avgF(1,ind+2)];
            
            [~, ind2] = min(abs(ymeasT(yk) - yankT));
            yapproxT = [yapproxT yankT(ind2)];
            ynu = [ynu yanknu(ind2)];
            yankid = [yankid ind2];
        end
    end
    
%% Produce a cell array of the fluorescence signal across each force range (for equilibrium constant calculation)

% This fills the cell array over the entire trace, not just the selected
% time interval. Uses the averaged fluorescence trace (from plotting).
for yy = 1:length(allyanks)-1
    meanT = mean(avgF(t_trap>allyanks(yy) & t_trap<allyanks(yy+1)));
    full_fl_sig = apdr1(t_apd>allyanks(yy) & t_apd<allyanks(yy+1));
    yr = {str2double(Date) datafnumber allyanks(yy) allyanks(yy+1) meanT oligo_conc full_fl_sig};
    YankResults = [YankResults; yr];
end

    
%% Set up a plot for analyzing binding/unbinding events

    if fluor_plots || step_plots
        
        figure('NumberTitle', 'off','FileName', [SaveFigDir Date '_' num2str(datafnumber, '%03d')...
            '_' probe], 'Name', [Date '_' num2str(datafnumber, '%03d') ' Analysis'])
        
        if fluor_plots
            if step_plots
                subplot(511)
            else
                subplot(411)
                box on
            end
            hold on
            
            apdrMax = max(apdav(t_apdav>plotinterval(1) & t_apdav<plotinterval(2)));
            apdrMin = min(apdav(t_apdav>plotinterval(1) & t_apdav<plotinterval(2)));
            t_apdin = t_apd(t_apd>plotinterval(1) & t_apd<plotinterval(2));
            t_apdinav = t_apdav(t_apdav>plotinterval(1) & t_apdav<plotinterval(2));
            apdrin = apdr1(t_apd>plotinterval(1) & t_apd<plotinterval(2));
            apdrinav = apdav(t_apdav>plotinterval(1) & t_apdav<plotinterval(2));

%             plot(t_apdin, apdrin, 'Color', Cy3)
            plot(t_apdinav, apdrinav, 'Color', Cy3)
            ylim([(apdrMin-apdrMin*0.1) (apdrMax+apdrMax*0.1)+.2])

            if length(t_yanking)>=2
                plot(t_yanking(2:end), mean(apdav,2)'+.1, '+m', 'MarkerSize', 10, 'linew', 2)
            end
            if ~isempty(t_unyanking)
                plot(t_unyanking, mean(apdav,2)'+.1, 'xm', 'MarkerSize', 10, 'linew', 2)
            end
            if ~isempty(t_binding)
                plot(t_binding, mean(apdav,2)'+.1, '+b', 'MarkerSize', 10, 'linew', 2)
            end
            if ~isempty(t_unbinding)
                plot(t_unbinding, mean(apdav,2)'+.1, 'xb', 'MarkerSize', 10, 'linew', 2)
            end
            ylabel('Photon Rate (kHz)')
            title([probe '  ' Date '\_' num2str(datafnumber)])

            if plotinterval(2)>t_trap(end)
                xlim([plotinterval(1) t_trap(end)])
            else
                xlim([plotinterval(1) plotinterval(2)])
            end
            ax2 = gca;
            set(gca,'XTickLabel', [])
        end
        
        if fluor_plots
            if step_plots
                subplot(512)
            else
                subplot(412)
                box on
            end
            hold on
            
%             apdravediffin = apdravediff(t_apd>plotinterval(1) & t_apd<plotinterval(2));
            apdravediffin = apdravediff;

            plot(t_apdin, apdravediffin)
%             plot(t_apd(1:end-floor(fwinsize*ApdBW)), -log10(fprob))
%             apdravediffMax = max(apdravediff(t_apd>plotinterval(1) & t_apd<plotinterval(2)));
%             apdravediffMin = min(apdravediff(t_apd>(plotinterval(1)+2) & t_apd<plotinterval(2)));
            apdravediffMax = max(apdravediff);
            apdravediffMin = min(apdravediff);
            % Guides
%             plot([t_apd(1) t_apd(end)], fpcutoff*[1 1], '--r')
            plot([t_apdin(1) t_apdin(end)], fltransthresh*[1 1], '--r')
            plot([t_apdin(1) t_apdin(end)], -fltransthresh*[1 1], '--r')
            
            ylim([(apdravediffMin+apdravediffMin*.1) (apdravediffMax+apdravediffMax*.1)])
            ylabel('Fluor averaged')
            ax2 = [ax2 gca];
            set(gca,'XTickLabel',[])
        end

        if fluor_plots
            if step_plots
                subplot(513)
            else
                subplot(413)
                box on
            end
            hold on

            plot(t_local, dextLocal-mean(dextLocal), 'k')
        elseif step_plots && ~fluor_plots
            subplot(311)
            hold on
            plot(t_local, dextLocal-mean(dextLocal), 'k')
            if ~isempty(t_binding)
                plot(t_binding, 2, '+k', 'MarkerSize', 10, 'linew', 2)
            end
            if ~isempty(t_unbinding)
                plot(t_unbinding, 2, 'xk', 'MarkerSize', 10, 'linew', 2)
            end
            title([probe '  ' Date '\_' num2str(datafnumber)])
        end
        ylabel('\DeltaExtension (nm)')
        ax2 = [ax2 gca];
        set(gca,'XTickLabel',[])
        
        if step_plots
            if fluor_plots
                subplot(514)
            else
                subplot(312)
                box on
            end
            hold on
            
            windowPts = floor(desired_trap_bw*windowSize);

            m = t_local(1:end-windowPts);
            plot(t_local(1:end-windowPts), -log10(prob), 'k')
            thresh = zeros(length(t_local),1);
            thresh(thresh==0) = probcutoff;
            plot(t_local, thresh, '--b')
            ylabel('t test results')
            ax2 = [ax2 gca];
            set(gca,'XTickLabel',[])
        end
        
        if fluor_plots
            if step_plots
                subplot(515)
            else
                subplot(414)
                box on
            end
        elseif step_plots && ~fluor_plots
            subplot(313)
        end
        hold on

        plot(t_local,AF_local(1,:),'b',t_local,-BF_local(1,:),'r',t_local,F_local,'k')
        ylabel('Force (pN)')
        xlabel('Time (s)')
        ax2 = [ax2 gca];
        
        linkaxes(ax2, 'x')
    end
    
%% Analyze binding events
    
    alldeltaextbind = [];
    for ie = 1:nbinding
        tevent = t_binding(ie);
        
        % Check what happened before event: see if enough time to measure
        % trap
        bindbefore = t_binding(t_binding < tevent);
        if isempty(bindbefore)
            bindbefore = plotinterval(1); % Changed 120402
        end
        closestbinding = min(tevent - bindbefore);
        
        unbindbefore = t_unbinding(t_unbinding <= tevent);
        if isempty(unbindbefore)
            unbindbefore = plotinterval(1); % Changed 150224
        end
        closestunbinding = min(tevent - unbindbefore);
        bg = mean(apdr1(ceil(unbindbefore(end)*ApdBW)+1:floor(tevent*ApdBW))); % [kHz] Background signal. Need +1 in case unbindbefore = 0.
        if length(DataSets{ii}) >= 12
            bgnorm = bg * 90/FBSP * 10^(NE/10)/10^2; % Background signal normalized to NE=20, green FB QPD setpoint=90 mV.
        else
            bgnorm = 0;
        end
        
        yankbefore = t_yanking(t_yanking <= tevent);
        if isempty(yankbefore)
            yankbefore = 0;
        end
        [closestyanking, yind] = min(tevent - yankbefore);
        stopBeforeYankBefore = closestyanking - 0.5; % Yank times are in the middle of the transition. Need to cut analysis before yank. -KW 110722
        
        unyankbefore = t_unyanking(t_unyanking <= tevent);
        if isempty(unyankbefore)
            unyankbefore = 0;
        end
        closestunyanking = min(tevent - unyankbefore);
        stopBeforeUnyankBefore = closestunyanking - 0.5; % Unyank times are in the middle of the transition. Need to cut analysis before unyank. -KW 110722
        
        mintbeforebinding = min([closestbinding; closestunbinding; closestyanking; closestunyanking]);
        availptbefore = mintbeforebinding * bw_trap - 1; % Need to move 1 point away, since this point might be affected by signal averaging.
%         availptbefore = mintbeforebinding * bw_trap-2;
        
        % Note that recording these as rate constants is exactly the WRONG
        % way to do it. 1/<tau> is NOT the same as <1/tau>. NO NO NO NO NO.
        % Should be recording LIFETIMES, not RATES. Thankfully, these are
        % immediately re-inverted once you get to plot_cut_results3, so
        % everything works out anyway. 161004
        if mintbeforebinding == closestunbinding
            lifeOn = closestunbinding; % [s]
            if strcmp(probe,'12mer_good') && str2double(Date)>=120725
                conc = -1; % I don't know what the concentration of these were.
            else
                conc = oligo_conc; % [nM]
            end
            rateOn = 1/lifeOn; % [1/s]
            kon = rateOn/conc; % [1/(nM*s)]
%             kon = rateOn/bgnorm; % [1/(kHz*s)]
            tau_ub = lifeOn * conc; % [nM s]
        else
            conc = oligo_conc;
            kon = -1; % Negative numbers will be cut later in plot_cut_results3.
            tau_ub = -1;
        end
        
%         useptsbefore = floor(min([trapdpoints; availptbefore]));
        useptsbefore = floor(min([avgptcutoff; availptbefore]));

        % Check what happened after event: see if enough time to measure trap
        bindafter = t_binding(t_binding > tevent);
        if isempty(bindafter)
            bindafter = 1e6;
        end
        closestbindingafter = min(bindafter - tevent);
        
        unbindafter = t_unbinding(t_unbinding > tevent);
        if isempty(unbindafter)
            unbindafter = 1e6;
        end
        closestunbindingafter = min(unbindafter - tevent);
        
        yankafter = t_yanking(t_yanking > tevent);
        if isempty(yankafter)
            yankafter = 1e6;
        end
        closestyankingafter = min(yankafter - tevent);
        stopBeforeYankAfter = closestyankingafter - 0.4; % Yank times are in the middle of the transition. Need to cut analysis before yank begins. -KW 110722
        
        unyankafter = t_unyanking(t_unyanking > tevent);
        if isempty(unyankafter)
            unyankafter = 1e6;
        end
        closestunyankingafter = min(unyankafter - tevent);
        stopBeforeUnyankAfter = closestunyankingafter - 0.4; % Unyank times are in the middle of the transition. Need to cut analysis before unyank begins. -KW 110722
        
        closetoend = min(min([plotinterval(2) t_trap(end)]) - tevent);
        
        mintafterbinding = min([closestbindingafter; closestunbindingafter; stopBeforeYankAfter; stopBeforeUnyankAfter; closetoend]);
%         availptafter = mintafterbinding * bw_trap - 2;
        availptafter = mintafterbinding * bw_trap - 1; % Need to move 1 point away, since this point might be affected by signal averaging.
        
%         useptsafter = floor(min([trapdpoints; availptafter]));
        useptsafter = floor(min([avgptcutoff; availptafter]));
        
        % For events that are too short, don't get an extension change. -KW 120124
        if ~(availptbefore > avgptcutoff && availptafter > avgptcutoff)
            continue % This event is no good, skip to next
        end
        
        t_before = t_trap(t_trap < tevent);
        beforeindices = ((length(t_before)-useptsbefore):(length(t_before))-1);
        t_before1 = t_before(beforeindices);
        
        % Important: we need to add 1 here because the points to either
        % side of the transition can be affected by the downsampling. Also
        % remember that we found the events by fluorescence, which is now
        % at a different frequency.
        t_after = t_trap(t_trap > tevent);
        afterindices = (1:useptsafter) + 1;
        t_after1 = t_after(afterindices);
        
        % Added this section 120629. Was previously only in "unbinding"
        % section.
        if ~isempty(ynu) && removeyankoffsets % Remove yank offsets from position and force data.
            if ~FB
                dext_before = dext(yankid(yind)+1, t_trap <= tevent);
                dext_after = dext(yankid(yind)+1, t_trap > tevent);
            else
                dext_before = dext(1, t_trap <= tevent);
                dext_after = dext(1, t_trap > tevent);
            end
            f_Abefore = Aforce(yankid(yind)+1, beforeindices);
            f_Bbefore = Bforce(yankid(yind)+1, beforeindices);
            
            f_Aafter0 = Aforce(yankid(yind)+1, t_trap > tevent);
            f_Bafter0 = Bforce(yankid(yind)+1, t_trap > tevent);
        else % No yank
            dext_before = dext(1, t_trap <= tevent);
            f_Abefore = Aforce(1, beforeindices);
            f_Bbefore = Bforce(1, beforeindices);
            
            dext_after = dext(1, t_trap > tevent);
            f_Aafter0 = Aforce(1, t_trap > tevent);
            f_Bafter0 = Bforce(1, t_trap > tevent);
        end
        
        f_Aafter = f_Aafter0(afterindices);
        f_Bafter = f_Bafter0(afterindices);
        meanforce = (mean(abs(f_Abefore))+mean(abs(f_Bbefore))+mean(abs(f_Aafter))+mean(abs(f_Bafter)))/4;
        diffF = mean([f_Aafter f_Abefore]) + mean([f_Bafter f_Bbefore]);
        
        %             dext_before = dext(1, t_trap <= tevent);
        dext_before1 = dext_before(beforeindices);
        extmeanbefore = mean(dext_before1);
        rmsbefore = rms(dext_before1 - extmeanbefore);
        
        %             dext_after = dext(1, t_trap > tevent);
        dext_after1 = dext_after(afterindices);
        extmeanafter = mean(dext_after1);
        rmsafter = rms(dext_after1 - extmeanafter);
        deltaextbind = extmeanafter - extmeanbefore;

        % If looking at labeled/unlabeled binding events, check if this
        % event is in the list for "unlabeled" events.
        if step_find && length(DataSets{ii})>=9
            stepind = find(t_all==t_binding(ie));
            isunlabel = find(unlabeled==stepind);
            if isempty(isunlabel)
                label = 1;
            else
                label = 0;
            end
        else
            label = 1;
        end
        if OT % No labeled binding events on the old trap.
            label = 0;
        end
        
        alldeltaextbind = [alldeltaextbind; [str2double(Date) datafnumber ie meanforce deltaextbind rmsbefore rmsafter mintbeforebinding mintafterbinding cal.kappaA cal.kappaB cal.alphaA cal.alphaB tevent mean([f_Aafter f_Abefore]) mean([f_Bafter f_Bbefore]) diffF label kon bg bgnorm conc]]; %#ok<*AGROW>
        
        % Plot steps found by t test
        if fluor_plots && step_plots
            subplot(511)
            if label
                plot(tevent, mean(apdr1,2)'+1, '+r', 'MarkerSize', 10, 'linew', 2)
            else
                plot(tevent, mean(apdr1,2)'+1, '+k', 'MarkerSize', 10, 'linew', 2)
            end
            subplot(513)
        elseif ~fluor_plots && step_plots
            subplot(311)
            plot(tevent, 2, '+k', 'MarkerSize', 10, 'linew', 2)
            subplot(311)
        end
        
        if fluor_plots && ~step_plots
            subplot(413)
        end
        if fluor_plots || step_plots
            plot(t_before1, dext_before1-mean(dextLocal), 'r+', t_after1, dext_after1-mean(dextLocal), 'b+')
        end
    end
    
    AllResults.binding = [AllResults.binding; alldeltaextbind];
    
%% Analyze unbinding events

    alldeltaextunbind = [];
    lifetimes = [];
    ruptures = [];
    for ie = 1:nunbinding
        tevent = t_unbinding(ie);
        % Determine what happened most recently before the unbinding event:
        % Another unbinding event? Can't measure lifetime.
        % A yank? Can't measure lifetime.
        % A binding event? Lifetime from binding, no trap adjust,
        % determine if there was a yank before unbinding.
        
        % Check what happened before event: See if enough time to measure
        % trap, get lifetime.
        bindbefore = t_binding(t_binding <= tevent);
        if isempty(bindbefore)
            bindbefore = plotinterval(1);
        end
        closestbinding = min(tevent - bindbefore);
        sig = mean(apdr1(ceil(bindbefore(end)*ApdBW)+1:floor(tevent*ApdBW))); % [kHz] Signal size
        
        Nph = 0;
        if strcmp(probe,'fluorLife')
            bgfl = mean(apdr1(ceil(tevent*ApdBW)+30:ceil(tevent*ApdBW)+50)); % Background fluorescence, averaged up to 50 points after signal drop.
            bgSubSig = apdr1(ceil(bindbefore(end)*ApdBW):floor(tevent*ApdBW)) - bgfl;
            Nph = sum(bgSubSig.*1000)/ApdBW;
            sig = sig - mean(apdr1(ceil(tevent*ApdBW):ceil(tevent*ApdBW)+50)); % Subtract some background signal
            
            if fluor_plots
                subplot(411)
                hold on
                plot(tevent+30/apdavBW, mean(apdav,2)'+.1, '+k', 'MarkerSize', 10, 'linew', 2)
                plot(tevent+50/apdavBW, mean(apdav,2)'+.1, 'xk', 'MarkerSize', 10, 'linew', 2)
            end
        end
        
        unbindbefore = t_unbinding(t_unbinding < tevent);
        if isempty(unbindbefore)
            unbindbefore = plotinterval(1);
        end
        closestunbinding = min(tevent - unbindbefore);
        
        yankbefore = t_yanking(t_yanking <= tevent);
        if isempty(yankbefore)
            yankbefore = plotinterval(1);
        end
        [closestyanking, yind] = min(tevent - yankbefore);
        stopBeforeYankBefore = closestyanking - 0.5; % Yank times are in the middle of the transition. Need to cut analysis before yank. -KW 110722
        if strcmp(probe,'12mer0mg_rupt') || strcmp(probe,'9mer_rupt')
            stopBeforeYankBefore = closestbinding;
        end
        
        unyankbefore = t_unyanking(t_unyanking <= tevent);
        if isempty(unyankbefore)
            unyankbefore = plotinterval(1);
        end
        closestunyanking = min(tevent - unyankbefore);
        stopBeforeUnyankBefore = closestunyanking - 0.5; % Unyank times are in the middle of the transition. Need to cut analysis before unyank. -KW 110722
        
        mintbeforebinding = min([closestbinding; closestunbinding; stopBeforeYankBefore; stopBeforeUnyankBefore]);
%         availptbefore = mintbeforebinding * bw_trap - 2;
        availptbefore = mintbeforebinding * bw_trap - 1; % Need to move 1 point away, since this point might be affected by signal averaging.
%         useptsbefore = floor(min([trapdpoints; availptbefore]));
        useptsbefore = floor(min([avgptcutoff; availptbefore]));
                
        % Check what happened after event: see if enough time to measure trap.
        bindafter = t_binding(t_binding > tevent);
        if isempty(bindafter)
            bindafter = plotinterval(2);
        end
        closestbindingafter = min(bindafter - tevent);
        
        unbindafter = t_unbinding(t_unbinding > tevent);
        if isempty(unbindafter)
            unbindafter = 1e6;
        end
        closestunbindingafter = min(unbindafter - tevent);
        
        yankafter = t_yanking(t_yanking > tevent);
        if isempty(yankafter)
            yankafter = 1e6;
        end
        closestyankingafter = min(yankafter - tevent);
        stopBeforeYankAfter = closestyankingafter - 0.5; % Yank times are in the middle of the transition. Need to cut analysis before yank begins. -KW 110722
        
        unyankafter = t_unyanking(t_unyanking > tevent);
        if isempty(unyankafter)
            unyankafter = 1e6;
        end
        closestunyankingafter = min(unyankafter - tevent);
        stopBeforeUnyankAfter = closestunyankingafter - 0.5; % Unyank times are in the middle of the transition. Need to cut analysis before unyank begins. -KW 110722
        
        mintafterbinding = min([closestbindingafter; closestunbindingafter; stopBeforeYankAfter; stopBeforeUnyankAfter]);
%         availptafter = mintafterbinding * bw_trap - 2;
        availptafter = mintafterbinding * bw_trap - 1; % Need to move 1 point away, since this point might be affected by signal averaging.
%         useptsafter = floor(min([trapdpoints; availptafter]));
        useptsafter = floor(min([avgptcutoff; availptafter]));  
                
        t_before = t_trap(t_trap < tevent);
        beforeindices = ((length(t_before)-useptsbefore):(length(t_before))-1);
        t_before1 = t_before(beforeindices);
        
        t_after = t_trap(t_trap > tevent);
        if useptsafter > 0
            afterindices = (1:useptsafter) + 1;
            afterfluorind = (1:useptsafter*2) + 1;
        else
            afterindices = 1;
        end
        if afterindices(end) > length(t_after) % This handles case at end of data
            ir = afterindices(end) - length(t_after);
            afterindices = afterindices(1:end-ir);
        end
        t_after1 = t_after(afterindices);
        
        if ~isempty(ynu) && removeyankoffsets % Remove yank offsets from position and force data.
            if ~FB
                dext_before = dext(yankid(yind)+1, t_trap <= tevent);
                dext_after = dext(yankid(yind)+1, t_trap > tevent);
            else
                dext_before = dext(1, t_trap <= tevent);
                dext_after = dext(1, t_trap > tevent);
            end
            f_Abefore = Aforce(yankid(yind)+1, beforeindices);
            f_Bbefore = Bforce(yankid(yind)+1, beforeindices);
            
            f_Aafter0 = Aforce(yankid(yind)+1, t_trap > tevent);
            f_Bafter0 = Bforce(yankid(yind)+1, t_trap > tevent);
        else % No yank
            dext_before = dext(1, t_trap <= tevent);
            f_Abefore = Aforce(1, beforeindices);
            f_Bbefore = Bforce(1, beforeindices);
            
            dext_after = dext(1, t_trap > tevent);
            f_Aafter0 = Aforce(1, t_trap > tevent);
            f_Bafter0 = Bforce(1, t_trap > tevent);
        end

        f_Aafter = f_Aafter0(afterindices);
        f_Bafter = f_Bafter0(afterindices);
        %         meanforce = (mean(abs(f_Abefore))+mean(abs(f_Bbefore))+mean(abs(f_Aafter))+mean(abs(f_Bafter)))/4;
        %         meanforcebefore = (mean(abs(f_Abefore))+mean(abs(f_Bbefore)))/2;
        %         meanforceafter = (mean(abs(f_Aafter))+mean(abs(f_Bafter)))/2;
        meanforce = (mean(f_Abefore)-mean(f_Bbefore)+mean(f_Aafter)-mean(f_Bafter))/4; % Changed 140207
        meanforcebefore = (mean(f_Abefore)-mean(f_Bbefore))/2;
        meanforceafter = (mean(f_Aafter)-mean(f_Bafter))/2;
        if meanforce < 0
            disp('Mean force went negative!')
        end
        
        % Get rupture force
        thresh = 0.1; % [s]
        ruptF_A = Aforce(1, t_trap<tevent+thresh & t_trap>tevent-thresh);
        ruptF_B = Bforce(1, t_trap<tevent+thresh & t_trap>tevent-thresh);
        ruptF = (ruptF_A(1) - ruptF_B(1))/2;
        
        % Get loading rate when rupture happens
        if ~isempty(f_Abefore) && ~isempty(f_Bbefore)
            F_beforerupt = (f_Abefore(end) - f_Bbefore(end))/2;
            F_afterrupt = (f_Aafter(1) - f_Bafter(1))/2;
            loadrate = (F_afterrupt - F_beforerupt)/(t_after1(1) - t_before1(end)); % [pN/s]
            tafterind = find(t_trap == t_after1(1));
            pullrate = -(traponly(tafterind) - traponly(beforeindices(end)))*XAODCal/(t_after1(1) - t_before1(end)); % [nm/s]
        else
            F_beforerupt = 0;
            F_afterrupt = 0;
            loadrate = 0;
            tafterind = 0;
            pullrate = 0;
        end

        dext_before1 = dext_before(beforeindices);
        extmeanbefore = mean(dext_before1);
        rmsbefore = rms(dext_before1 - extmeanbefore);
        
        dext_after1 = dext_after(afterindices);
        extmeanafter = mean(dext_after1);
        rmsafter = rms(dext_after1 - extmeanafter);
  
        deltaextunbind = extmeanafter - extmeanbefore;
        
        yesyank = closestyanking < closestbinding || closestunyanking < closestbinding; % Then a yank happened during the event
        if closestyankingafter < 1 || closestunyankingafter < 1 % If unbinding event is too close to a yank, can't be sure whether fluor drop was from unbinding or result of yank. -KW 120221
            life = -1;
        else
            if yesyank && closestbinding-1>=plotinterval(1) && closestbinding+1<=plotinterval(1)
                if closestyanking < closestunyanking
                    if closestbinding-closestyanking>4 && closestbinding-closestyanking<5 % Remove events which have too long or short of a gap between binding and yanking/unyanking.
                        life = closestyanking;
                        lifeforce = meanforceafter;
                        lifetimes = [lifetimes; [str2double(Date) datafnumber ie lifeforce life cal.kappaA cal.kappaB cal.alphaA cal.alphaB tevent deltaextunbind -1 -1 -1 -1 label]];
                    else
                        life = -1;
                        lifeforce = meanforceafter;
                    end
                else
                    if closestbinding-closestunyanking>4 && closestbinding-closestunyanking<5
                        life = closestunyanking;
                        lifeforce = meanforceafter;
                        lifetimes = [lifetimes; [str2double(Date) datafnumber ie lifeforce life cal.kappaA cal.kappaB cal.alphaA cal.alphaB tevent deltaextunbind -1 -1 -1 -1 label]];
                    else
                        life = -1;
                        lifeforce = meanforceafter;
                    end
                end
            elseif ~yesyank && bindbefore(end)~=plotinterval(1) && bindbefore(end)>unbindbefore(end)% Fixed. -KW 120221
                life = closestbinding;
                lifeforce = meanforceafter;
                lifetimes = [lifetimes; [str2double(Date) datafnumber ie lifeforce life cal.kappaA cal.kappaB cal.alphaA cal.alphaB tevent deltaextunbind sig Nph NE FBSP label]];
            else % If data starts with probe bound, don't record any lifetime. -KW 110705.
                life = -1;
                lifeforce = meanforceafter;
            end
        end
        ruptures = [ruptures; ruptF loadrate pullrate];
        
        % Plot lifetime
        str = sprintf('%3.1f s', life);
        if mod(ie,2)
            textoff = 1;
        else
            textoff = -1;
        end
        if fluor_plots && step_plots
            subplot(511)
%             text(t_unbinding(ie), mean(apdr1(beforeindices), 2)' + textoff, str);
        elseif fluor_plots && ~step_plots
            subplot(411)
            if strcmp(probe,'12mer0mg_rupt') || strcmp(probe,'9mer_rupt')
                text(t_unbinding(ie), mean(apdr1(beforeindices),2)'+textoff, ['F_c = ' num2str(ruptF,'%0.2f') ' pN'])
            end
%             text(t_unbinding(ie), mean(apdr1(beforeindices), 2)' + textoff, str);
            subplot(413)
        elseif ~fluor_plots && step_plots
            subplot(311)
            text(t_unbinding(ie), 3, str);
        end
        
        % If looking at labeled/unlabeled binding events, check if this
        % event is in the list for "unlabeled" events.
        if step_find && length(DataSets{ii})>=9
            stepind = find(t_all==t_unbinding(ie));
            isunlabel = find(unlabeled==stepind);
            if isempty(isunlabel)
                label = 1;
            else
                label = 0;
            end
        elseif ~step_find
            label = 1;
        end
        
        if availptbefore > avgptcutoff && availptafter > avgptcutoff % Trap is OK
            
            alldeltaextunbind = [alldeltaextunbind; [str2double(Date) datafnumber ie meanforce deltaextunbind rmsbefore rmsafter life useptsbefore useptsafter cal.kappaA cal.kappaB cal.alphaA cal.alphaB tevent label]];
            
        % Plot steps found by t test
        if fluor_plots && step_plots
            subplot(511)
            if label
                plot(tevent, mean(apdr1,2)'+1, 'xr', 'MarkerSize', 10, 'linew', 2)
            else
                plot(tevent, mean(apdr1,2)'+1, 'xk', 'MarkerSize', 10, 'linew', 2)
            end
            subplot(513)
        elseif ~fluor_plots && step_plots
            subplot(311)
            plot(tevent, 2, 'xk', 'MarkerSize', 10, 'linew', 2)
        end

        if fluor_plots || step_plots
            plot(t_before1, dext_before1-mean(dextLocal), 'g+', t_after1, dext_after1-mean(dextLocal), 'm+')
        end
        end
    end
    clear dext Aforce Bforce % Save memory
    
    AllResults.unbinding = [AllResults.unbinding; alldeltaextunbind];
    AllResults.lifetimes = [AllResults.lifetimes; lifetimes];
    AllResults.rupture = [AllResults.rupture; ruptures];

%     evalin('base', 'clear(calfilename)')
%     evalin('base', 'clear(offsetfilename)') % Clear out stuff to save memory
end

% Save some RAM by clearing out temp variables
clear alldeltaextbind alldeltaextunbind
toc