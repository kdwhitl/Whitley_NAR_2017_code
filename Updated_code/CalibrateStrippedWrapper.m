% Author: Matt Comstock
% Date last modified: 150728 (added argument 'startpath' in case this is
% used for the old trap).

% This function returns two structure variables by assigning them to the
% command window (the 'base'). These are "allcal" and "allfit", which are
% saved as "cal[Date]_[filenum]" and "fit[Date]_[filenum]".

function CalibrateStrippedWrapper(Date,rootfile,showplots,writecal,startpath)
if nargin == 0
    Date = '161202';
    datafnumber = 33;
    rootfile = [Date '_' num2str(datafnumber,'%03d') '.dat'];
    
    prompts = {'beadA (nm)','beadB (nm)','fsamp (Hz)','fXYHigh (Hz)', ...
        'fSumHigh (Hz)','fLow (Hz)','Averaging windows'};
    defaults = {'880','900','125000','62500','2000','100','100'};
    param = getnumbers('Enter parameters:',prompts,defaults);
end

if nargin >= 4
else
    showplots = 1;
    writecal = 1;
end

global datadirectories
global laptop
global beaddiameters

if nargin < 5
    startpath = [datadirectories Date '\'];
end

calfilename = ['cal' rootfile(1:(end-4))];
fitfilename = ['fit' rootfile(1:(end-4))];

data = ReadMattFile_Wrapper(startpath,rootfile);

%which signal sets to analyze
calsets(1) = 1;%Trap A
calsets(2) = 1;%Trap B
calsets(3) = 0;%Detection laser

%Strep beads are usually 790 nm and ADig (protein G) beads are usually 860 nm

param = zeros(7,1);
param(7) = 200;%ave windows
param(5) = 8e3;%f sum high

if calsets(1) == 1 %get cal for trap 1 (A)
    param(1) = beaddiameters(1);
    % param(1) = 2100;%Bead A
    % param(1) = 790;%Bead A
    param(4) = 8e3;%fxy high
    param(6) = 600;%f low
    param(6) = 1.1e3;%f low. 170713
    fitmodel = 'AliasedHydro';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.A_X;
    data_sub.Y = data.A_Y;
    data_sub.Sum = data.A_Sum;
    
    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    allcal.beadA = param(1);
    allcal.alphaAX = cal.alphaAX;
    allcal.alphaAY = cal.alphaAY;
    allcal.kappaAX = cal.kappaAX;
    allcal.kappaAY = cal.kappaAY;
    allcal.AXoffset = cal.AXoffset;
    allcal.AYoffset = cal.AYoffset;
    
    allfit.fA = fit.f;
    allfit.AXSpecRaw = fit.AXSpecRaw;
    allfit.AYSpecRaw = fit.AYSpecRaw;
    allfit.predictedAX = fit.predictedAX;
    allfit.predictedAY = fit.predictedAY;
end

if calsets(2) == 1 %get cal for trap 2 (B)
    param(1) = beaddiameters(2);
    % param(1) = 2100;%Bead B
    % param(1) = 880;%Bead B
    param(4) = 8e3;%fxy high
    param(6) = 600;%f low
    param(6) = 1.1e3;%f low. 170713
    fitmodel = 'AliasedHydro';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.B_X;
    data_sub.Y = data.B_Y;
    data_sub.Sum = data.B_Sum;
    
    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    % All "cal" fields are named as if they were bead A from the
    % "CalibrateExact_Stripped" function, but they need to be renamed here,
    % since this is actually the calibration for bead B.
    allcal.beadB = param(1);
    allcal.alphaBX = cal.alphaAX;
    allcal.alphaBY = cal.alphaAY;
    allcal.kappaBX = cal.kappaAX;
    allcal.kappaBY = cal.kappaAY;
    allcal.BXoffset = cal.AXoffset;
    allcal.BYoffset = cal.AYoffset;
    
    allfit.fB = fit.f;
    allfit.BXSpecRaw = fit.AXSpecRaw;
    allfit.BYSpecRaw = fit.AYSpecRaw;
    allfit.predictedBX = fit.predictedAX;
    allfit.predictedBY = fit.predictedAY;
end

if calsets(3) == 1 %get cal for detection laser on trap
    param(1) = 1000;%Bead B
    param(4) = 15e3;%fxy high
    param(6) = 200;%f low
    fitmodel = 'AliasedFiltered';
    %fitmodel = 'AliasedHydro';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.C_X;
    data_sub.Y = data.C_Y;
    data_sub.Sum = data.C_Y_Sum;
    
    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    allcal.beadC = param(1);
    allcal.alphaCX = cal.alphaAX;
    allcal.alphaCY = cal.alphaAY;
    allcal.kappaCX = cal.kappaAX;
    allcal.kappaCY = cal.kappaAY;
    allcal.CXoffset = cal.AXoffset;
    allcal.CYoffset = cal.AYoffset;
    
    allfit.fC = fit.f;
    allfit.CXSpecRaw = fit.AXSpecRaw;
    allfit.CYSpecRaw = fit.AYSpecRaw;
    allfit.predictedCX = fit.predictedAX;
    allfit.predictedCY = fit.predictedAY;
    
end

assignin('base',calfilename,allcal);
assignin('base',fitfilename,allfit);

if showplots
    
    columns = sum(calsets);
    
    if laptop == 1
        if sum(calsets) == 3
            figure('Position',[222          89        1626         871])
        else
            figure('Position',[25          49        900         600])
        end
    else
        %figure('Position',[10   319   1400*columns/3   637])
    end
    
    set(gcf,'Name',['Calibration for ' rootfile(1:end-4)]);
    prev = 0;
    if calsets(1) == 1
        subplot(2,columns,1)
        loglog(allfit.fA, allfit.AXSpecRaw); hold on;
        loglog(allfit.fA, allfit.predictedAX, 'k');
        title(['AX: \kappa = ' num2str(allcal.kappaAX,3) ', \alpha = ' ...
            num2str(allcal.alphaAX,4) ', Offset = ' num2str(allcal.AXoffset,3)]);
        % ylim([1e-9 5e-8])
        ylabel('Power (V^2 s)')
        subplot(2,columns,1+columns)
        loglog(allfit.fA, allfit.AYSpecRaw); hold on;
        loglog(allfit.fA, allfit.predictedAY, 'k');
        title(['AY: \kappa = ' num2str(allcal.kappaAY,3) ', \alpha = ' ...
            num2str(allcal.alphaAY,4) ', Offset = ' num2str(allcal.AYoffset,3)]);
        % ylim([1e-9 5e-8])
        prev = 1;
        xlabel('Frequency (Hz)')
        ylabel('Power (V^2 s)')
    end
    
    if calsets(2) == 1
        subplot(2,columns,prev+1)
        loglog(allfit.fB, allfit.BXSpecRaw); hold on;
        loglog(allfit.fB, allfit.predictedBX, 'k');
        title(['BX: \kappa = ' num2str(allcal.kappaBX,3) ', \alpha = ' ...
            num2str(allcal.alphaBX,4) ', Offset = ' num2str(allcal.BXoffset,3)]);
        % ylim([1e-9 5e-8])
        subplot(2,columns,prev+1+columns)
        loglog(allfit.fB, allfit.BYSpecRaw); hold on;
        loglog(allfit.fB, allfit.predictedBY, 'k');
        % ylim([1e-9 5e-8])
        title(['BY: \kappa = ' num2str(allcal.kappaBY,3) ', \alpha = ' ...
            num2str(allcal.alphaBY,4) ', Offset = ' num2str(allcal.BYoffset,3)]);
        prev = prev + 1;
        xlabel('Frequency (Hz)')
    end
    
    if calsets(3) == 1
        subplot(2,columns,prev+1)
        loglog(allfit.fC, allfit.CXSpecRaw); hold on;
        loglog(allfit.fC, allfit.predictedCX, 'k');
        title(['CX: \kappa = ' num2str(allcal.kappaCX,3) ', \alpha = ' ...
            num2str(allcal.alphaCX,4) ', Offset = ' num2str(allcal.CXoffset,3)]);
        subplot(2,columns,prev+1+columns)
        loglog(allfit.fC, allfit.CYSpecRaw); hold on;
        loglog(allfit.fC, allfit.predictedCY, 'k');
        title(['CY: \kappa = ' num2str(allcal.kappaCY,3) '. \alpha = ' ...
            num2str(allcal.alphaCY,4) ', Offset = ' num2str(allcal.CYoffset,3)]);
    end
    
end

if (isempty(laptop) || laptop ~= 1) && writecal == 1 %means running on instrument computer
    if showplots
        saveas(gcf,[startpath rootfile(1:(end-4)) '_cal.fig']);
    end
    %attempting to get the data to labview automatically
    calout = [allcal.alphaAX allcal.alphaBX allcal.kappaAX allcal.kappaBX allcal.AXoffset allcal.BXoffset];
    fid2 = fopen([startpath 'current.cal'],'w','ieee-be');
    fwrite(fid2,calout,'float64');
    fclose(fid2);
end