% Author: Matt Comstock
% Date created: ?
% Date last modified: 140908 by Kevin Whitley

% Plot results of hybrid step size analysis
SaveFigDir = 'C:\Users\Kevin\Documents\UHR_Fleezer_Data\Oligo_data'; % on laptop
if isdir([SaveFigDir '\' probe])
    SaveFigDir = [SaveFigDir '\' probe];
end

todaysDate = datestr(now,'yymmdd');
FigSize = [25   49  900     600];

%% Data sets

plot_ext = 0;
plot_life = 1;
plot_onrate = 0;
plot_T = 0;
plot_rupture = 0;
plot_Nphotons = 0;
plot_bg_conc = 0;
plot_normSig = 0;
plot_Kd = 0;
plot_ton_conc = 0;

date_cut = 0; % Cut on date when the data was taken.

label_only = 1;
unlabel_only = 0;
removePB = 0; % Remove photobleaching events by looking at unbinding extension change.
pbcutoff = 0; % [nm] Remove any lifetimes whose corresponding extension change is smaller than this.
removeShort = 0; % Remove extension changes from short events.
shortcutoff = 0; % [s] Remove any extension changes whose corresponding lifetimes are smaller than this.
cut_fluor_back = 0;
pb_correct = 0; % Include photobleaching rate for off rates
sig_cut = 0; % Cut on signal size
sigthresh = 400; % Hz
cut_tether_noise = 0;
extmat = [];

extbins = -3:0.2:3; % Bins for extension plots.
lowbound = 2; % [s] Lowest bound for fitting lifetime exponentials.
on_lowbound = 20; % [s nM]

lives = [];
meanT = [];
errort_ub = [];
tau1 = [];
errortau1 = [];
avph = [];
stph = [];
F_yank = [];
Kd = [];
Nphsave = {};

for Tset = 1:8
switch Tset
    case 11
        force = 1;
    case 1
        force = 3;
    case 2
        force = 5;
    case 3
        force = 8;
    case 4
        force = 10;
    case 5
        force = 13;
    case 6
        force = 15;
    case 7
        force = 17.5;
    case 8
        force = 20;
    case 9
        force = 22.5;
    case 10
        force = 25;
end
cut = [force-1 force+1];
if plot_T
    cut = [0 25]; % all
end
SaveName = [SaveFigDir '\' todaysDate '_' probe '_' num2str(force,'%2.0f') 'pN_'];

%% Cut data

% Cut on force range
Data_Binding =          cutmat(AllResults.binding, 4, cut);
Data_Unbinding =        cutmat(AllResults.unbinding, 4, cut);
if ~isempty(AllResults.lifetimes)
    Data_Lifetimes =    cutmat(AllResults.lifetimes, 4, cut);
end

if plot_Kd && ~isempty(YankResults)  % This is a cell array, so can't use cutmat.
    pp = 1;
    ymat = cell2mat(YankResults(:,1:6));
    
    [new_ymat, indy] = unique(ymat,'rows'); % Remove redundant intervals/traces
    yforce = new_ymat(:,5);
    forceind = find(yforce>cut(1) & yforce<cut(2));
    F_yank(Tset) = mean(yforce(yforce>cut(1) & yforce<cut(2)));
    
    YRmat = YankResults(indy,:); % Take out redundant traces
    Data_Yank = YRmat(forceind,:); % Only this force interval
end

% Cut on fluorescence background
% if cut_fluor_back
%     bkcut = [0 10]; % [kHz]
%     Data_Binding = cutmat(Data_Binding, 20, bkcut);
% 
%     figure('Name',[probe '_' force 'pN_background_fluorescence'])
%     hist(Data_Binding(:,20))
%     xlabel('Background signal (Hz)')
%     ylabel('Counts')
% end

% Cut on difference between A and B forces
if ~isempty(Data_Binding) && 0
    Data_Binding =      cutmat(Data_Binding, 17, [0 4]);
end

% For ext, cut on tether noise
if cut_tether_noise
    bnoisecut = [0 2];
    unbnoisecut = [0 1];
    Data_Binding =          cutmat(Data_Binding, 6, bnoisecut);
    Data_Binding =          cutmat(Data_Binding, 7, bnoisecut);
    Data_Unbinding =        cutmat(Data_Unbinding, 6, unbnoisecut);
    Data_Unbinding =        cutmat(Data_Unbinding, 7, unbnoisecut);
end

% For ext, cut on trap stiffness (really just important for 9merCy3 data)
bkappaAcut = [0.1 0.6];
ubkappaAcut = [0.1 0.8];
%         bkappaAcut = [0.175 0.225]; % For old 9merCy3 data at 3pN.
%         ubkappaAcut = [0.175 0.225];
Data_Binding =          cutmat(Data_Binding,10,bkappaAcut);
Data_Unbinding =        cutmat(Data_Unbinding,11,ubkappaAcut);

% Cut on time of event in time trace
if 0
    teventcut = [0 500];
    Data_Binding =      cutmat(Data_Binding,14,teventcut);
    Data_Unbinding =    cutmat(Data_Unbinding,15,teventcut);
    Data_Lifetimes =    cutmat(Data_Lifetimes,10,teventcut);
end

% Remove photobleaching events for high forces (>8pN)
if removePB && ~isempty(Data_Lifetimes)
    if Tset~=1 && Tset~=11
        Data_Lifetimes = cutmat(Data_Lifetimes, 11, [pbcutoff max(Data_Lifetimes(:,11))]);
    else
        Data_Lifetimes = cutmat(Data_Lifetimes, 11, [min(Data_Lifetimes(:,11)) pbcutoff]);
    end
end

% Cut on date
if date_cut
%     datecut = [110101 150701];
%     datecut = [150601 160821];
    datecut = [160401 160801];
%     datecut = [120320 150904];
%     datecut = [130905 130908];
    Data_Lifetimes =    cutmat(Data_Lifetimes,1,datecut);
    Data_Binding =      cutmat(Data_Binding,1,datecut);
    Data_Unbinding =    cutmat(Data_Unbinding,1,datecut);
    
    % Cut on individual time trace
%     tracecut = [143 143];
%     Data_Lifetimes = cutmat(Data_Lifetimes,2,tracecut);
%     Data_Binding = cutmat(Data_Binding,2,tracecut);
%     Data_Unbinding = cutmat(Data_Unbinding,2,tracecut);
end

% Cut on fluorescence background
if cut_fluor_back
    bkcut = [0 20]; % [kHz]
    Data_Binding =      cutmat(Data_Binding, 20, bkcut);

    figure('Name',[probe '_' force 'pN_background_fluorescence'])
    hist(Data_Binding(:,20))
    xlabel('Background signal (Hz)')
    ylabel('Counts')
end
        
% Cut labeled/unlabeled events (for half-and-half t test)
if label_only
    labelrange = [1 1];
elseif unlabel_only
    labelrange = [0 0];
else
    labelrange = [0 1];
end
if ~isempty(Data_Binding) && length(Data_Binding(1,:))>=18 && (~isempty(Data_Lifetimes) && length(Data_Lifetimes(1,:))>=16)
    Data_Binding =          cutmat(Data_Binding, 18, labelrange);
    Data_Unbinding =        cutmat(Data_Unbinding, 16, labelrange);
    Data_Lifetimes =        cutmat(Data_Lifetimes, 16, labelrange);
end

if ~isempty(Data_Lifetimes)
    lives = Data_Lifetimes(:,5);
end

% Remove extension changes from short events
if removeShort && ~isempty(Data_Lifetimes)
    Data_Lifetimes = [Data_Lifetimes lives];
    Data_Lifetimes =    cutmat(Data_Lifetimes, 5, [shortcutoff max(Data_Lifetimes(:,5))]);
    Data_Lifetimes(:,[5,11]) = Data_Lifetimes(:,[11,5]);
    Data_Unbinding = Data_Lifetimes;
end

if removePB && removeShort
    disp('Do not cut photobleaching and short events at the same time!')
    t = pouiyoui576s; % Will (appropriately) give an error and stop parsing.
end

%% A few more things before plotting

% Add in events manually (couldn't find them automatically for some reason,
% so they were never added into AllResults).
others = [];
otherF = [];
if strcmp(probe,'12mer1') || strcmp(probe,'12mer_good') % Pretty sure we need this 'good' part in there.
    if  Tset==5
        others = [209.4];
        otherF = [12.13];
    elseif Tset==6
        others = [92.7 68.7 575 7 329 97 97 121 185 352 977 703.5 157 612.5 294.5 461.5];
        otherF = [15.04 14.79 14.4 14.4 15.3 15.3 15.3 15.3 14.9 14.35 14.35 15.91 14.58 14.58 14.85 14.85];
    elseif Tset==7
        others = [200 125.3 174.4 53 22.4 183.3 91 285 102 32 138.5 20 156 50 237];
        otherF = [17.48 17.56 17.87 17.87 17.41 17.98 17.5 17.5 18.09 18.09 18.09 18.09 18.09 18.09 18.09];
    elseif Tset==8
        others = [82.9 55.6];
        otherF = [20.29 19.61];
    end
end
if strcmp(probe,'10mer') || strcmp(probe,'10mer_cleanset')
    switch Tset
        case 4
            others = [163.8 99.8 60.2 31 49.9 70 41 22];
            otherF = [10.16 10.16 10.16 10.16 10.04 10.07 10.07 10.05];
        case 1
            others = [135 193 92 29 56.7 58];
            otherF = [2.9 2.9 3.16 3.16 3.08 3.08];
        case 2
            others = [201 56 138];
            otherF = [5.08 5.08 5.18];
        case 3
            others = [108.5 73.8 119.4 27 47 30.4 40.2];
            otherF = [8.01 8.01 8.02 8.02 8.02 7.94 7.94];
    end
end
if strcmp(probe,'rna10mer')
    switch Tset
        case 2
            others = 164.5;
            otherF = 5.07;
        case 3
            others = 152;
            otherF = 8.19;
        case 4
            others = [40 123.1];
            otherF = [9.75 10.15];
        case 5
            others = [25 72.8 1.5 73.1 15.9 59.7 21.5 17.3];
            otherF = [11.72 12.37 12.38 12.38 12.38 12.38 12.43 12.43];
    end
end
% if strcmp(probe,'rna11mer') % With gloxy
%     if Tset==5
%         others = 193;
%         otherF = 12.3;
%     end
% end
if strcmp(probe,'seq3')
    if Tset==4
        others = 172;
        otherF = 10;
    elseif Tset==5
        others = 42;
        otherF = 12;
    end
end
if strcmp(probe,'rna8mer')
    % Removed all this due to new beads, 130822
%     if Tset==1
%         others = [5.7 22.7 11 73.4 26.6 8.1]';
%         otherF = [3 3 3 3 3 3]';
%     elseif Tset==2
%         others = [11 6.5 1.8 3.4 5.2]';
%         otherF = [5 5 5 5 5]';
    if Tset==3
        others = 17.2;
        otherF = 8;
%         others = 4.5;
%         otherF = 8;
    end
end
if strcmp(probe,'12mer0mg')
    if Tset==5
        others = [161.9 53.5 36 69.9];
        otherF = [12.85 12.54 12.56 12.78];
    end
    if Tset==6
        other = [112.8 4.5 32.6 50.5];
        otherF = [15.14 15.14 15.14 15];
    end
    if Tset==7
        other = 30.3;
        otherF = 17.88;
    end
end
if strcmp(probe,'9mer0mg')
    if Tset==1
        others = [42 50];
        otherF = [3.3 3];
    end
    if Tset==11
        others = 18;
        otherF = 1.5;
    end
end
if ~isempty(lives)
    lives = [lives; others'];
end
if strcmp(probe,'10mer0mg')
    if Tset==11
        others = [185.8 202.4 176.8 31];
        otherF = [1.5 1.5 1.5 1.5];
    end
    if Tset==1
        others = 162;
        otherF = 3;
    end
    if Tset==2
        others = 106;
        otherF = 5;
    end
    if Tset==3
        others = 50;
        otherF = 8;
    end
end
forces = [Data_Binding(:,4); otherF'];

display(' ');
meanT(Tset) = mean(forces);
errorT = std(forces) / sqrt(length(forces));
display(['Mean binding tension = ' num2str(meanT(Tset), '%0.2f') ' +/- ' num2str(errorT,'%0.2f') ' pN'])

%% Plot delta extension - transparent plot
if (~isempty(Data_Binding) || ~isempty(Data_Unbinding)) && plot_ext
    
    meanunbinding = [];
    figure('NumberTitle', 'off', 'Name', [probe '  ' num2str(meanT(Tset),'%0.1f')...
        'pN Oligo Hybridization Extension Changes'], 'Pos', FigSize,'FileName', [SaveName 'dext'])
    hold on
    legendstr = {};
    
    disp(['N_b = ' num2str(length(Data_Binding(:,5)),'%02d') ', N_ub = ' num2str(length(Data_Unbinding(:,5)),'%02d')])
    if ~isempty(Data_Binding)
        hist(Data_Binding(:,5), extbins)
        hb = hist(Data_Binding(:,5), extbins);
        h = get(gca, 'chi');
        set(h(1), 'facea', 0.5, 'facecolor', 'b')
    
        meanbinding = mean(Data_Binding(:,5));
        errorbinding = std(Data_Binding(:,5)) / sqrt(length(Data_Binding(:,5)));
        display(['Mean binding D_Ext = ' num2str(meanbinding, '%0.2f') ' +/- ' num2str(errorbinding,'%0.1g') ' nm'])

        legendstr = {'Binding'};
    else
        hb = [];
    end
    
    if ~isempty(Data_Unbinding)
        hist(Data_Unbinding(:,5), extbins)
        hunb = hist(Data_Unbinding(:,5), extbins);
        h = get(gca, 'chi');
        set(h(1), 'facea', 0.5, 'facecolor', 'g')

        meanunbinding = mean(Data_Unbinding(:,5));
        errorunbinding = std(Data_Unbinding(:,5)) / sqrt(length(Data_Unbinding(:,5)));
        display(['Mean unbinding D_ext = ' num2str(meanunbinding, '%0.2f') ' +/- ' num2str(errorunbinding,'%0.1g') ' nm'])
        
        meanboth = mean([-Data_Unbinding(:,5); Data_Binding(:,5)]);
        errboth = std([-Data_Unbinding(:,5); Data_Binding(:,5)]) / sqrt(length(Data_Unbinding(:,5))+length(Data_Binding(:,5)));
        display(['Mean both Dext = ' num2str(meanboth, '%0.2f') ' +/- ' num2str(errboth,'%0.1g') ' nm'])
        
        legendstr{length(legendstr) + 1} = 'Unbinding';
    else
        hunb = [];
    end

    xlabel('\Deltax_h (nm)')
    ylabel('Counts')
    title([probe '. ' num2str(meanT(Tset),'%0.2f') ' pN'])
    if ~isempty(meanunbinding)
        text((extbins(end)-extbins(1))/20+extbins(1),max([hb hunb])/1.3,{['\Delta_b = ' num2str(meanbinding,'%0.2f') ' +/- ' num2str(errorbinding,'%0.1g') ' nm'] ['\Delta_{unb} = ' num2str(meanunbinding,'%0.2f') ' +/- ' num2str(errorunbinding,'%0.1g') ' nm']})
    else
        text((extbins(end)-extbins(1))/20+extbins(1),max([hb hunb])/1.3,{['\Delta_b = ' num2str(meanbinding,'%0.2f') ' +/- ' num2str(errorbinding,'%0.1g') ' nm']})
    end
    legend(legendstr, 'Location', 'Northwest', 'Color', 'none')
    xlim([extbins(1) extbins(end)])
    
    allexts = [Data_Binding(:,5); -Data_Unbinding(:,5)];
    extmat = [extmat; allexts ones(length(allexts),1)*meanT(Tset) ones(length(allexts),1)*Tset];
%     c1 = length(allexts);
%     [r1, c2] = size(extmat);
%     if isempty(extmat)
%         extmat = [extmat allexts];
%     elseif c1 > r1;
%         extmatnan = ones(c1,c2).*NaN;
%         extmatnan(1:r1,1:c2) = extmat;
%         extmat = [extmatnan allexts];
%     else
%         allexts = [allexts; ones(r1-c1,1).*NaN];
%         extmat = [extmat allexts];
%     end
    
%     figure
%     hist([Data_Binding(:,5); -Data_Unbinding(:,5)], extbins)
%     hunb = hist([Data_Binding(:,5); -Data_Unbinding(:,5)], extbins);
%     h = get(gca, 'chi');
%     set(h(1), 'facea', 0.5, 'facecolor', 'r')
%     xlabel('\Deltax_h (nm)')
%     ylabel('Counts')
%     title([probe '. ' num2str(meanT(Tset),'%0.2f') ' pN'])
%     xlim([extbins(1) extbins(end)])
    
    % Plot noise distribution
    if 0
        figure('NumberTitle','off','Name',[forceCase 'pN_noise'],'Pos',[12    62   600   600])
        subplot(211)
        hist([Data_Binding(:,6);Data_Binding(:,7)])
        xlabel('Sigma (nm)')
        ylabel('Counts')
        title('Binding \DeltaExt Noise Dist.')
        subplot(212)
        hist([Data_Unbinding(:,6);Data_Unbinding(:,7)])
        xlabel('Sigma (nm)')
        ylabel('Counts')
        title('Unbinding \DeltaExt Noise Dist.')
    end
end

%% Plot oligo lifetimes
if ~isempty(Data_Lifetimes) && plot_life
    
    % Get signal sizes for events, subtracting their backgrounds
    Data_Binding = sortrows(Data_Binding);
    Data_Lifetimes = sortrows(Data_Lifetimes);
    
    if length(Data_Lifetimes(1,:)) >=12
        sig = Data_Lifetimes(:,12) .* 1000; % [Hz] Fluorophore signals + background
        
        % Subtract the background signal for each file individually (excitation
        % settings may have varied from one to the next).
        [nfisig, is] = unique(Data_Lifetimes(:,1:2),'rows'); % Data files for lifetimes
        is = [is; length(Data_Lifetimes(:,1)) + 1]; % Add in last elements
        
        % Need to match up files in Data_Lifetimes and Data_Binding (one has
        % background data, other has signal data).
        filecut = [];
        for pp = 1:size(nfisig)
            filecutind = find(nfisig(pp,1)==Data_Binding(:,1) & nfisig(pp,2)==Data_Binding(:,2));
            filecut = [filecut; Data_Binding(filecutind,:)];
        end
        [nfiles, ia] = unique(filecut(:,1:2),'rows'); % Data files involved
        ia = [ia; length(filecut(:,1)) + 1]; % Add in last elements
        
        % Need to make sure ia and is are the same size (all lifetimes have
        % backgrounds to match).
        [~, indf] = ismember(nfisig, nfiles,'rows');
        is(indf==0) = [];
        
        sigtot = [];
        bk = [];
        for nn = 1:length(ia)-1
            bk(nn) = mean(filecut(ia(nn):ia(nn+1)-1,20)) .* 1000; % [Hz] Background signal in each file.
            sigtot = [sigtot; sig(is(nn):is(nn+1)-1) - bk(nn)];
        end
        if sig_cut % Cut on signal size
            lives = Data_Lifetimes(:,5); % No signal size info on manually-added lifetimes, so remove.
            sigtot = sigtot(lives>=lowbound);
            lives = lives(sigtot<sigthresh);
            sigtot = sigtot(sigtot<sigthresh);
        end
        sigall = mean(sigtot);
    end
    
    %     figure
    %     plot(1:length(lives),lives,'.k')
    lives = lives(lives>=lowbound);
    tau1(Tset) = mean(lives) - lowbound;
    errortau1(Tset) = std(lives) / sqrt(length(lives));
    display(['Lifetime data points = ' num2str(length(lives))])
    display(['Mean lifetime = ' num2str(tau1(Tset), '%0.2f') ' +/- ' num2str(errortau1(Tset), '%0.1g') ' s'])
    
    iqb = iqr(lives); % Interquartile range for exponential distribution
    binsize = 2*iqr(lives)*length(lives)^(-1/3); % Freedman-Diaconis rule for determining optimal bin sizes.
    
    figure('NumberTitle', 'off', 'Name', [probe '  ' num2str(meanT(Tset),'%0.1f')...
        'pN Oligo Hybridization Lifetimes'], 'Pos', FigSize, 'FileName', [SaveName 'lifetimes'])
    hold on
    box on
    xlim([0 max(lives)+binsize])
    xlabel('t_b (s)')
    ylabel('Survival probability')
    title([probe ' lifetimes at ' num2str(meanT(Tset),'%0.2f') ' pN'])
    set(gca,'YTick', 0:0.2:1)
    
    cbinsize = 0.001;
%     cbins = 0:cbinsize:max(lives);
    cbins = 0:cbinsize:180; % for plotting
    
    li = histc(lives,cbins);
    cdist = cumsum(li)/max(cumsum(li));
    cbinsfit = cbins(cbins > lowbound);
    cdistfit = cdist(cbins > lowbound);
    
    cdistfit2 = cdistfit(diff(cdistfit)~=0);
    cbins2 = cbinsfit(diff(cdistfit)~=0);
    %     lifesort = sort(lives);
    %     cdist2 = cumsum(lifesort);
    
    %     plot(cbins2', 1-cdistfit2, 'b', 'linew', 2)
    plot(cbinsfit', 1-cdistfit, 'm', 'linew', 2)

    if pb_correct % Correct for photobleaching
        kph = 2.22e-5; % [1/photons]
        expmodel = @(koff,b) exp(-(b-lowbound).*(koff + kph.*sigall));
        [efit2, Rd, Jd] = nlinfit(cbins2', 1-cdistfit2, expmodel, 1./tau1(Tset));
    else
        expmodel = @(koff,b) exp(-(b-lowbound).*koff);
        [efit2, R2, J2] = nlinfit(cbins2', 1-cdistfit2, expmodel, 1./tau1(Tset));
    end
    
    [efit, R, J] = nlinfit(cbinsfit', 1-cdistfit, expmodel, 1./tau1(Tset));
    [~, se] = nlparci2(efit, R, 'Jacobian', J);
    expect = expmodel(efit,cbins2);
    residsq = sum((1-cdistfit2' - expect).^2./expect);
    %     [h, p] = chi2gof(1-cdistfit' - expect);
    
    %     modelBins = lowbound:0.1:max(lives)+binsize;
    modelBins = lowbound:0.1:180; % for plotting
%     plot(modelBins, expmodel(efit,modelBins), '--k', 'linew', 2)
    plot(modelBins, expmodel(efit2,modelBins), '--r', 'linew', 2)
    text((cbins(end)-cbins(1))/2,0.5, {['\tau_b = ' num2str(1/efit2,'%0.2f') ' +/- ' num2str(1/efit2/sqrt(length(lives)),'%0.2g') ' s']})
%     display(['Fit life = ' num2str(efit,'%0.2f') ' +/- ' num2str(efit/sqrt(length(lives)),'%0.2f') ' s'])
    display(['Fit life2 = ' num2str(1/efit2,'%0.2f') ' +/- ' num2str(1/efit2/sqrt(length(lives)), '%0.1g') ' s'])
%     display(['Residual square = ' num2str(residsq,'%0.2f')])
    
    % Plot lifetime histogram (with bins)
    if 0
        % First bin is just empty. Want binning to start with low bound. Since
        % matlab chooses all bin widths based on the size of the first bin,
        % need to make that one the same size as all the rest, even if it goes
        % into negative values.
        figure
        bins = [lowbound-binsize lowbound:binsize:max(lives)+binsize];
        hh = histc(lives, bins);
        hh = hh(1:end - 1);
        xvals = (bins(1:end-1) + (bins(2)-bins(1))/2)';
        bar(xvals, hh, 'b')
        
%         Remove zeros (or not)
        binsfitrange = [lowbound max(lives)];
        hhfit = hh(hh > 0 & xvals > binsfitrange(1) & xvals < binsfitrange(2));
        binsfit = xvals(hh > 0 & xvals > binsfitrange(1) & xvals < binsfitrange(2));
        hhfit = hh(xvals > binsfitrange(1) & xvals < binsfitrange(2));
        binsfit = xvals(xvals > binsfitrange(1) & xvals < binsfitrange(2));
        
%         Fit to weighted single exponential using 'fit'
        [expfit, gof] = fit(binsfit, hhfit, 'exp1', 'Weights', hhfit, 'Lower', [1 -Inf], 'Upper', [1000 0], 'StartPoint', [length(lives) -1/tau1(Tset)]);
        plot(modelBins, expfit.a * exp(modelBins*expfit.b), '--c', 'linew', 2)
        display(['Fit life = ' num2str(-1/expfit.b,'%0.2f') ' s'])
    end
    
    % Fit and plot cumulative distribution
    if 0
        % Fit to cumulative distribution function (single exponential)
        cumexp = @(mu, x) 1-exp(-(x-lowbound)./mu);
        [mufit, R, J] = nlinfit(cbinsfit', cdistfit./max(cdist), cumexp, tau1(Tset));
        [~, muerr] = nlparci2(mufit,R,'Jacobian',J);
        display(['Cum mu = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.1g') ' s'])
        %     plot(modelBins, sum(hh).*exp(-modelBins/mufit), 'k--', 'linew', 2)
        %     text((bins(end)-bins(1))/2,max(hh)/2-5, {['\tau_{cf} = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.2f') ' s']});
        
        figure
        hold on
        plot(cbins,cdist./max(cdist),'k')
        plot(modelBins, cumexp(mufit, modelBins),'--r','linew',2)
        text((cbins(end)-cbins(1))/2, 1/2, {['\tau_{cf} = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.1g') ' s']})

        xlabel('Lifetime (s)')
        ylabel('Cumulative probability')
    end
end

%% Plot 'unbound' lifetimes
if ~isempty(Data_Binding) && plot_onrate
    
    % These were (incorrectly) recorded as 'on' rate constants, so need to
    % invert. Remember that false measurements were all recorded as -1, so
    % need to take >0.
    t_ub = 1./Data_Binding(Data_Binding(:,19)>0, 19);
%     t_ub = Data_Binding(Data_Binding(:,19)>0, 19);
    
    onlife = t_ub(t_ub >= on_lowbound);
    tau_ub(Tset) = mean(onlife) - on_lowbound; % [nM s]
    errort_ub(Tset) = std(t_ub) / sqrt(length(t_ub));
    
    display(['N_ub = ' num2str(length(t_ub), '%02d')])
    display(['<t_ub> = ' num2str(tau_ub(Tset), '%2.2e') ' +/- ' num2str(errort_ub(Tset), '%0.1g') ' nM s'])

    figure('NumberTitle', 'off', 'Name', [probe '  ' num2str(force,'%2.0f')...
        'pN Oligo unbound state lifetimes'], 'FileName', [SaveName, 't_ub'])
    hold on
    box on
    
    iqbon = iqr(t_ub); % Interquartile range for exponential distribution
    binsize = 2*iqr(t_ub)*length(t_ub)^(-1/3); % Freedman-Diaconis rule for determining optimal bin sizes.
    
    cbinsize = 1;
    cbins = 0:cbinsize:max(onlife);
    
    li = histc(onlife,cbins);
    cdist = cumsum(li)/max(cumsum(li));
    cbinsfit = cbins(cbins > on_lowbound);
    cdistfit = cdist(cbins > on_lowbound);
    
    cdistfit2 = cdistfit(diff(cdistfit)~=0);
    cbins2 = cbinsfit(diff(cdistfit)~=0);
    
%     plot(cbins2', 1-cdistfit2, 'b', 'linew', 2)
    plot(cbinsfit', 1-cdistfit, 'm', 'linew', 2)

    expmodel = @(tauf,b) exp(-(b-on_lowbound)./tauf);
%     efit = nlinfit(cbinsfit', 1-cdistfit, expmodel, tau_ub(Tset));
    efit2 = nlinfit(cbins2', 1-cdistfit2, expmodel, tau_ub(Tset));
    
    modelBins = on_lowbound:1:max(onlife)+binsize;
%     plot(modelBins, expmodel(efit,modelBins), '--k', 'linew', 2)
    plot(modelBins, expmodel(efit2,modelBins), '--r', 'linew', 2)
    text((cbins(end)-cbins(1))/2,0.5, {['\tau_{ub} = ' num2str(efit2,'%0.2g') ' +/- ' num2str(efit2/sqrt(length(onlife)),'%0.1g') ' s']})
    display(['Fit on life = ' num2str(efit2,'%0.2f') ' +/- ' num2str(efit2/sqrt(length(onlife)),'%0.1g') ' s'])
    
    if 0
        figure
        ubbins = 0:binsize:(max(t_ub)+binsize);
        
        hh = histc(t_ub, ubbins);
        hh = hh(1:end - 1);
        xvals = (ubbins(1:end-1) + (ubbins(2)-ubbins(1))/2)';
        bar(xvals, hh, 'b')
        xlim([0 ubbins(end)])
        
        title([probe ' unbound time at ' num2str(force,'%2.1f') ' pN'])
        
        % Remove zeros
        ubbinsfit = [5e-3 1e4]; % First bin is a bit off due to finite length of time traces.
        hhfit = hh(xvals > ubbinsfit(1) & xvals < ubbinsfit(2));
        binsfit = xvals(xvals > ubbinsfit(1) & xvals < ubbinsfit(2));
        
        % Fit to weighted single exponential using 'fit'
        [expfit, gof] = fit(binsfit, hhfit, 'exp1', 'Weights', hhfit, 'Lower', [1 -Inf], 'Upper', [1000 0], 'StartPoint', [10 -1/tau_ub(Tset)]);
        
        errfit = -1/(expfit.b*sqrt(sum(hhfit)));
        display(['Fit ub life = ' num2str(-1/expfit.b,'%2.2e') ' +/- ' num2str(-1/expfit.b/sqrt(length(t_ub)),'%0.1g') ' nM s'])
        
        hold on
        modelBins = ubbins(1):1:ubbins(end);
        plot(modelBins, expfit.a * exp(modelBins*expfit.b), '--c', 'linew', 2)
        %     plot(modelBins, singlexp(taufitw,modelBins), 'k--', 'linew', 2)
        
        text((ubbins(end)-ubbins(1))/2,max(hh)/2, {['\tau_{ub} = ' num2str(tau_ub(Tset),'%2.2e') ' +/- ' num2str(errort_ub(Tset),'%0.1g') ' nM s']})% ['\tau_f = ' num2str(lifefit,'%0.2f') ' +/- ' num2str(errorfit, '%0.2f') ' s']});
    end

%     li = histc(t_ub(t_ub>on_lowbound),cbins);
%     cdist = cumsum(li);
%     cbinsfit = cbins(cbins > on_lowbound);
%     cdistfit = cdist(cbins > on_lowbound);
%     cumexp = @(mu, x) 1-exp(-(x-on_lowbound)./mu);
%     [mufit, R, J] = nlinfit(cbinsfit', cdistfit./max(cdist), cumexp, tau_ub(Tset));
%     [~, muerr] = nlparci2(mufit,R,'Jacobian',J);
%     display(['Cum mu ub = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.2f') ' nM s'])
%     plot(modelBins, sum(hh).*exp(-modelBins/mufit), 'k--', 'linew', 2)
%     text((bins(end)-bins(1))/2,max(hh)/2-5, {['\tau_{cf} = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.2f') ' s']});
    
    % Plot cumulative distribution
    if 0
        figure
        hold on
        plot(cbins,cdist./max(cdist),'k')
        plot(modelBins, cumexp(mufit, modelBins),'--r','linew',2)
        text((cbins(end)-cbins(1))/2, 1/2, {['\tau_{ub, cf} = ' num2str(mufit,'%0.2f') ' +/- ' num2str(muerr,'%0.1g') ' s']})
        ylim([0 1])

        xlabel('Lifetime (s)')
        ylabel('Cumulative probability')
    end
end

%% Plot dext vs tether time (scatter)
if 0
    figure('NumberTitle', 'off', 'Name', [forceCase 'PN_ext_vs_time'], 'Pos', [20 100 600 600])
    plot(Data_Binding(:,end), Data_Binding(:,5), 'b.')
end

%% Plot lifetimes vs tether time
if 0
    figure('NumberTitle', 'off', 'Name', [forceCase 'pN_lifetime_vs_time'], 'Pos', [650 100 600 600],...
        'FileName',[SaveFigDir '\' todaysDate '_' forceCase 'pN_lifetime_vs_time'])
    plot(Data_Lifetimes(:,end), Data_Lifetimes(:,5), 'b.')
    xlabel('Tether time (s)')
    ylabel('Lifetime (s)')
    
%     tevent = Data_Lifetimes(:,end);
%     timeBin = [];
%     
%     for times = 0:100:1000
%         avgLives = mean(lives(tevent > times & tevent < times+70));
%         timeBin = [timeBin; avgLives, times];
%     end
%     figure
%     plot(timeBin(:,2)+35, timeBin(:,1), '*r', 'MarkerSize', 10);
end

%% Plot calibration values (mostly just important for 9merCy3 data)
if 0
    figure('Pos', [25    62   841   619])
    subplot(2, 1, 1)
    hist(Data_Binding(:,10), 0:0.005:0.25)
    xlabel('\kappa_{AX} (pN/nm)')
    ylabel('Counts')
    title('Binding Calibration Dist.')
    
    subplot(2, 1, 2)
    hist(Data_Unbinding(:,11), 0:0.005:0.25)
    xlabel('\kappa_{BX} (pN/nm)')
    ylabel('Counts')
    title('Unbinding Calibration Dist.')
    
    % bkaxstats = vec_stats(Data_Binding(:,10));
    bkaxstats = mean(Data_Binding(:,10));
    bkaxerror = std(Data_Binding(:,10))/sqrt(length(Data_Binding(:,10)));
    display(['Mean binding k_ax = ' num2str(bkaxstats, '%0.3f') ' +/- ' num2str(bkaxerror,'%0.3f') ' pN/nm'])
    % bkbxstats = vec_stats(Data_Binding(:,11));
    bkbxstats = mean(Data_Binding(:,11));
    bkbxerror = std(Data_Binding(:,11))/sqrt(length(Data_Binding(:,11)));
    display(['Mean binding k_bx = ' num2str(bkbxstats, '%0.3f') ' +/- ' num2str(bkbxerror,'%0.3f') ' pN/nm'])
end

%% Plot correlation between step-finder and fluorescence
if 0
    [rowS, ~] = size(Data_stepBinding);
    [row, ~] = size(Data_Binding);
    finalDeltaBind = [];
    for ee = 1:rowS
        bindstepTime = Data_stepBinding(ee,end);
        deltaBind = [];
        for ff = 1:row
            bindTime = Data_Binding(ff,end);
            deltaBind(ff) = bindstepTime - bindTime;
        end
        finalDeltaBind(ee) = min(abs(deltaBind));
    end
    
    [rowUS, ~] = size(Data_stepUnbinding);
    [rowU, ~] = size(Data_Unbinding);
    finalDeltaUnbind = [];
    for uu = 1:rowUS
        unbindstepTime = Data_stepUnbinding(uu,end);
        deltaUnbind = [];
        for nn = 1:rowU
            unbindTime = Data_Unbinding(nn,end);
            deltaUnbind(nn) = unbindstepTime - unbindTime;
        end
        finaDeltaUnbind(uu) = min(abs(deltaUnbind));
    end
    
    figure('Name', ['Stepfinder_fluorescence_correlation_' todaysDate], 'NumberTitle', 'off')
    hist(finalDeltaBind)
end

%% Plot tether tension
if plot_T
    T = mean(Data_Lifetimes(:,4));
    errorT = std(Data_Lifetimes(:,4)) / sqrt(length(Data_Lifetimes(:,4)));
    display(['Mean lifetime tension = ' num2str(T, '%0.1f') ' +/- ' num2str(errorT, '%0.1f') ' pN'])
    
    Tbins = 0.2;
    Fbins = cut(1):Tbins:cut(2);
    
    figure('Name', 'Oligo Hybridizing: Tension Analysis')
    hist(Data_Lifetimes(:,4), Fbins);
    hh = hist(Data_Lifetimes(:,4), Fbins);
%     xlim([Fbins(1) Fbins(end)])
    
    ht = text(Fbins(1)+0.2, max(hh)/2,['T = ' num2str(T, '%0.1f') ' +/- ' num2str(errorT, '%0.1f') ' pN']);
    
    xlabel('Tension (pN)')
    ylabel('Counts')
    set(gcf, 'FileName', SaveName);
end

%% Plot F_AX and F_BX difference distribution
if ~isempty(AllResults.binding) && 0
    T = mean(Data_Lifetimes(:,4));
    
    F_AX = cutmat(AllResults.binding,15,[-50 50]);
    F_BX = cutmat(AllResults.binding,16,[-50 50]);
    F_AX = F_AX(:,15);
    F_BX = F_BX(:,16);
    
    figure('Name',[probe '_' todaysDate ' (F_AX+F_BX)'])

    bin1 = 0:0.5:10;
    hist(F_AX+F_BX, bin1, 'b');
    xlabel('F_AX+F_BX (pN)')
    ylabel('Events')
end

%% Plot rupture force distribution
if plot_rupture
%     allrupt = cutmat(AllResults.rupture, 3, [1 1000]);
    allrupt = cutmat(AllResults.rupture, 3, [390 410]);
    
    figure('Name',[probe '_' todaysDate '_rupturedist'])
    bins = 0:2:40;
    hist(allrupt(:,1),bins)
    [rupth, edge] = hist(allrupt(:,1),bins');
    err_rupt = 1./sqrt(rupth);
    binwidth = bins(2)-bins(1);
    prob_rupt = rupth./(sum(rupth)*binwidth);
    bins = bins(prob_rupt > 0);
    err_rupt = err_rupt(prob_rupt > 0);
    prob_rupt = prob_rupt(prob_rupt > 0);
    
%     bar(bins, prob_rupt)
    xlabel('Rupture force (pN)')
    ylabel('Probability')
    
    if 1
        figure('Name',[probe '_' todaysDate '_loadingrates'])
        plot(allrupt(:,1),allrupt(:,2),'.k')
        xlabel('Rupture force (pN)')
        ylabel('Loading rate (pN/s)')
    end
    
    life_rupt = [];
    F_out = [];
    Fdot = [];
    for F_ind = 1:length(bins)
        kap_bead = 1/(1/mean(Data_Unbinding(:,11)) + 1/mean(Data_Unbinding(:,12))); % [pN/nm]
%         kap_bead = 0.005;
        kaptot = (1/kap_bead + 0.34*3263/XWLCContour(bins(F_ind),50,1000,4.14,1) + 0.6*10/XWLCContour(bins(F_ind),1,1000,4.14,1))^-1;
        vel = mean(allrupt(:,3)); % [nm/s]
        
%         Fdot(F_ind) = vel * (1/kap + (2*B*cont*lp*(1+B*bins(F_ind)*lp)) / (3+5*B*bins(F_ind)*lp + 8*(B*bins(F_ind)*lp)^(5/2)))^-1;
        Fdot(F_ind) = vel * kaptot;
        life_rupt(F_ind) = sum(prob_rupt(F_ind:end))*binwidth/(Fdot(F_ind)*prob_rupt(F_ind));
        F_out(F_ind) = bins(F_ind);
    end
    
    kapav = (1/kap_bead + 0.34*3263/XWLCContour(mean(allrupt(:,1)),50,1000,4.14,1) + 0.6*10/XWLCContour(mean(allrupt(:,1)),1,1000,4.14,1))^-1;
    Favdot = vel * kapav;
    avgtau = sqrt(pi/2)*std(allrupt(:,1))/Favdot;
%     hold on
%     plot(F_out,Fdot,'r')
    
    figure('Name',[probe '_' todaysDate '_rupture'])
    hold on
    errorbar(F_out,life_rupt,life_rupt.*err_rupt,'.k','linew',2,'MarkerSize',20)
    plot(mean(allrupt(:,1)),avgtau,'+b','linew',2,'MarkerSize',12)
    set(gca,'YScale','log')
    xlabel('Force (pN)')
    ylabel('\tau_b (s)')
    
    x1 = F_out(life_rupt > 0 & life_rupt < 1e6);
    y1 = life_rupt(life_rupt > 0 & life_rupt < 1e6);
    e1 = y1.*err_rupt(life_rupt > 0 & life_rupt < 1e6);
    nu = 2/3;
    kT = 4.14; % [pN nm]
    eqn = @(a,f) 1./(a(1).*(1-nu.*f.*a(3)./a(2)/kT).^(1/nu-1).*exp(a(2).*(1-(1-nu.*f.*a(3)./a(2)/kT).^(1/nu))));
    eqnw = @(a,f) 1./e1.*eqn(a,f); % Weighted equation
    [fitvals,R,J,covb,mse] = nlinfit(x1, y1./e1, eqnw, [1/20 5 .5]');
    Fplot = 0:0.5:35;
    plot(Fplot,eqn(fitvals,Fplot),'--r','linew',2)
    
    display(['Dx = ' num2str(fitvals(3),'%0.2f') ' nm'])
    display(['DG = ' num2str(fitvals(2),'%0.2f') ' kT'])
    display(['k(0) = ' num2str(fitvals(1),'%0.2e') ' 1/s'])
end

if plot_normSig
    Data_Binding = sortrows(Data_Binding);
    Data_Lifetimes = sortrows(Data_Lifetimes);
    sig = Data_Lifetimes(:,12) .* 1000; % [Hz] Fluorophore signals + background
    NEall = Data_Lifetimes(:,14);
    SPall = Data_Lifetimes(:,15);
    
    % Subtract the background signal for each file individually (excitation
    % settings may have varied from one to the next).
    [nfisig, is] = unique(Data_Lifetimes(:,1:2),'rows'); % Data files for lifetimes
    is = [is; length(Data_Lifetimes(:,1)) + 1]; % Add in last elements
    
    % Need to match up files in Data_Lifetimes and Data_Binding (one has
    % background data, other has signal data).
    filecut = [];
    for pp = 1:length(nfisig)
        filecutind = find(nfisig(pp,1)==Data_Binding(:,1) & nfisig(pp,2)==Data_Binding(:,2));
        filecut = [filecut; Data_Binding(filecutind,:)];
    end
    [nfiles, ia] = unique(filecut(:,1:2),'rows'); % Data files involved
    ia = [ia; length(filecut(:,1)) + 1]; % Add in last elements
    
    % Need to make sure ia and is are the same size (all lifetimes have
    % backgrounds to match).
    [~, indf] = ismember(nfisig, nfiles,'rows');
    is(indf==0) = [];
    
    sigtot = [];
    bk = [];
    ND = [];
    FBSP = [];
    for nn = 1:length(ia)-1
        bk(nn) = mean(filecut(ia(nn):ia(nn+1)-1,20)) .* 1000; % [Hz] Background signal in each file.
        sigtot = [sigtot; sig(is(nn):is(nn+1)-1) - bk(nn)];
        ND = [ND; NEall(is(nn):is(nn+1)-1)];
        FBSP = [FBSP; SPall(is(nn):is(nn+1)-1)];
    end
    signorm = sigtot .* 90./FBSP .* 10.^(ND./10)./10.^2; % Background signal normalized to NE=20, green FB QPD setpoint=90 mV.
%     signorm = sigtot .* 90./FBSP .* 10.^2./10.^(ND./10); % Background signal normalized to NE=20, green FB QPD setpoint=90 mV.
    
    figure('Name',[probe ' at ' num2str(mean(forces),'%0.1f') 'pN, signal size'], 'FileName', [SaveName '_sigsize'])
    hold on
    box on
    bins = 100:50:2000;
    h2 = histcounts(signorm,bins);
    %         h2 = h2(1:end-1);
    xvals = (bins(1:end-1) + (bins(2)-bins(1))/2);
    bar(xvals, h2)
    title(['Signal size after background subtraction for ' probe ...
        ' data. Signal normalized to normal excitation settings: ND = 20, SP = 90 mV.'])
    xlabel('Photon rate (kHz)')
    ylabel('Counts')
end

if plot_Nphotons
    
    figure('Name', [probe '_' todaysDate '_Nphotons'])
    Data_Binding = sortrows(Data_Binding);
    Data_Lifetimes = sortrows(Data_Lifetimes);
    sig = Data_Lifetimes(:,12) .* 1000; % [Hz] Fluorophore signals + background
    li2 = Data_Lifetimes(:,5); % [s] Event lifetimes
    
    % Subtract the background signal for each file individually (excitation
    % settings may have varied from one to the next).
    [nfisig, is] = unique(Data_Lifetimes(:,1:2),'rows'); % Data files for lifetimes
    is = [is; length(Data_Lifetimes(:,1)) + 1]; % Add in last elements
    
    % Need to match up files in Data_Lifetimes and Data_Binding (one has
    % background data, other has signal data).
    filecut = [];
    for pp = 1:size(nfisig)
        filecutind = find(nfisig(pp,1)==Data_Binding(:,1) & nfisig(pp,2)==Data_Binding(:,2));
        filecut = [filecut; Data_Binding(filecutind,:)];
    end
    [nfiles, ia] = unique(filecut(:,1:2),'rows'); % Data files involved
    ia = [ia; length(filecut(:,1)) + 1]; % Add in last elements
    
    % Need to make sure ia and is are the same size (all lifetimes have
    % backgrounds to match).
    [~, indf] = ismember(nfisig, nfiles,'rows');
    is(indf==0) = [];
    
    sigtot = [];
    bk = [];
    li3 = [];
    for nn = 1:length(ia)-1
        bk(nn) = mean(filecut(ia(nn):ia(nn+1)-1,20)) .* 1000; % [Hz] Background signal in each file.
        sigtot = [sigtot; sig(is(nn):is(nn+1)-1) - bk(nn)];
        li3 = [li3; li2(is(nn):is(nn+1)-1)]; % Need to match up lifetimes.
    end
%     sigtot = sig - mean(bg); % [Hz] Fluorophore signals
    Nph = sigtot .* li3; % (mean) Number of photons emitted for each event

    if sig_cut % Cut on signal size
        Nph = Nph(sigtot<sigthresh);
    end
    Nphsave{Tset} = Nph;
    avph(Tset) = mean(Nph);
    stph(Tset) = std(Nph)/sqrt(length(Nph));
    
    iqbon = iqr(Nph); % Interquartile range for exponential distribution
    binsize = 2*iqr(Nph)*length(Nph)^(-1/3); % Freedman-Diaconis rule for determining optimal bin sizes.
    
    lowbound = 500;
    Nph = Nph(Nph>lowbound);
    bins = [lowbound-binsize lowbound:binsize:max(Nph)+binsize];
    hh = histc(Nph, bins);
    hh = hh(1:end - 1);
    xvals = (bins(1:end-1) + (bins(2)-bins(1))/2)';
    bar(xvals, hh, 'r')
    
    % Fit to weighted single exponential using 'fit'
    [expfit, gof] = fit(xvals, hh, 'exp1', 'Weights', hh, 'Lower', [1 -Inf], 'Upper', [1000 0], 'StartPoint', [30 -1/30000]);
    
    errfit = -1/(expfit.b*sqrt(sum(hh)));
    display(['Mean Nph = ' num2str(avph(Tset),'%2.2e') ' +/- ' num2str(stph(Tset),'%2.2e') ' photons'])
    display(['Fit Nph = ' num2str(-1/expfit.b,'%2.2e') ' +/- ' num2str(-1/expfit.b/sqrt(length(Nph)),'%2.2e') ' photons'])
    
    hold on
    modelBins = bins(1):binsize/10:bins(end);
%     plot(modelBins, expfit.a * exp(modelBins*expfit.b), '--k', 'linew', 2)
% 
%     title([probe ' photons emitted before bleaching at ' num2str(mean(forces),'%0.1f') ' pN'])
%     xlabel('Number of photons')
%     ylabel('Counts')
%     display(['Avg photons at ' num2str(mean(forces),'%0.0f') 'pN: ' num2str(avph(Tset),'%0.0f') ' +/- ' num2str(stph(Tset),'%0.0f')])
    
    if 1
        figure('Name',[probe ' at ' num2str(mean(forces),'%0.1f') 'pN, signal size'])
        bins = 100:200:2000;
        h2 = histcounts(sigtot,bins);
%         h2 = h2(1:end-1);
        xvals = (bins(1:end-1) + (bins(2)-bins(1))/2);
        bar(xvals, h2)
        xlabel('Photon rate (kHz)')
        ylabel('Counts')
    end
end

% Plot background signal vs. concentration
if plot_bg_conc
    figure('Name', [probe '_' todaysDate '_background_vs_conc'], 'FileName', [SaveName 'bg_conc'])
    box on
    hold on
    
    [nfiles, inds] = unique(Data_Binding(:,1:2),'rows'); % Individual data files. Don't want repeats.
    
    conc = Data_Binding(inds, 22); % [nM]
    bgnorm = Data_Binding(inds, 21); % Background normalized to NE=20, green FB QPD setpoint=90 mV.
    
    % Remove zeros
    conc = conc(bgnorm > 0);
    bgnorm = bgnorm(bgnorm > 0);
    
    scatter(conc, bgnorm, '*k')
    xlabel('Concentration (nM)')
    ylabel('Normalized background (kHz)')
    
end
if plot_Kd && ~isempty(Data_Yank)
    figure('Name', [probe '_' todaysDate '_fluorsig_hist'], 'FileName', [SaveName 'fluorsig_hist'])
    box on
    hold on
    
    allfluor = [];
    [r, ~] = size(Data_Yank);
    FluorIntFact = 100;
    for nn = 1:r
        intfluor = apd_integrate(Data_Yank{nn,7},FluorIntFact) ./ FluorIntFact; % Not quite integrated. Would need bw for this.
        allfluor = [allfluor intfluor];
    end
    
    nbins = 50;
    binsize = (max(allfluor)-min(allfluor))/nbins;
    bins = min(allfluor):binsize:max(allfluor);
    [cts, edg] = histcounts(allfluor, bins);
    xvals = edg(1:end-1) + (edg(2)-edg(1))/2;
    bar(xvals, cts)
    
    midline = (max(allfluor)-min(allfluor))/2 + min(allfluor);
    plot(ones(1,max(cts)).*midline, 1:max(cts), '--r')
    xlabel('Photon rate (kHz)')
    ylabel('Counts')
    title([probe ' fluorescence signal at ' num2str(F_yank(Tset),'%0.1f') ' pN'])
    
    b_conc = length(allfluor(allfluor>midline)); % Amount of bound state
    unb_conc = length(allfluor(allfluor<midline)); % Amount of unbound state
    
    Kd(Tset) = unb_conc / b_conc;
    
    text(midline, max(cts), ['[unbound]/[bound] = ' num2str(Kd(Tset),'%0.2f')])
end

% Plot unbound-state lifetimes as a function of concentration. 160928
if plot_ton_conc
    figure('Name', [probe '_' todaysDate '_ton_vs_conc'], 'FileName', [SaveName 'ton_vs_conc'])
    box on
    hold on
    
    tauon = 1./Data_Binding(:,19); % [s nM]
    conc = Data_Binding(:,22); % [nM]
    
    ton = tauon./conc; % [s]
    conc = conc(ton>0);
    ton = ton(ton>0); % remove non-events
    
    if ~isempty(ton)
        plot(conc, ton, '+k')
        
        [cmean, tmean, ~, terr] = BoxcarAvg(conc, ton, 5);
        errorbar(cmean, tmean, terr, 'or', 'linew', 2, 'MarkerSize', 10)
        
        if length(cmean)>2
            newterr = terr./tmean.^2; % error propagation, 1/x (stdv)
            [m, R, J] = nlinfit(cmean, tmean, @(a,x)a(1)./x + a(2), [1 0], 'Weights', 1./newterr.^2);
            [~, se] = nlparci2(m, R, 'Jacobian', J);
            tc = 0:35;
            plot(tc, m(1)./tc+m(2), ':b')
%             text(5, 150, ['\tau_{on} = ' num2str(m(1)./1e9,'%2.1e') ' +/- ' num2str(se(1)./1e9,'%1.0e') ' s M'])
            text(5, 150, ['k_{on} = ' num2str(1e9/m(1),'%2.1e') ' +/- ' num2str(1e9*se(1)/m(1)^2,'%1.0e') ' s^{-1} M^{-1}'])
        end
    end
    xlabel('[oligo] (nM)')
    ylabel('t_{on} (s)')
    title([probe ' t_{on} vs. [oligo] at ' num2str(meanT(Tset),'%0.1f') ' pN'])

end
end
if plot_Nphotons
    % Convert cell array to matrix, filling in blanks with NaNs
    ms = max(cellfun('size',Nphsave,1));
    
    nphmat = [];
    for ii = 1:length(Nphsave)
        newcol = [Nphsave{ii}; ones((ms-length(Nphsave{ii})),1).*NaN];
        if any(newcol)
            nphmat = [nphmat newcol];
        end
    end
end
