% Author: Matt Comstock
% Date last modified: never (at least not by Kevin)

% This function appears to have come direcly from CalibrateExact7 (or some
% variant thereof). Only one bead calibration is done in this version,
% however, since the two beads are handled separately in
% CalibrateStrippedWrapped.

function [calibration,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,file,data)
%basic calibration routine that just generates the data, nothing else, for
%a single set of data
%090514 mjc

% Set Data Paths
global analysisPath;%folder where Analysis and Power folders are located

% fitmodel = 'AliasedHydro'; 
% fitmodel = 'AliasedFiltered';
% Set Parameter defaults
% if nargin < 2 %~exist('fsamp') | (exist('fsamp') & isempty(fsamp))
%     prompts = {'beadA (nm)','beadB (nm)','fsamp (Hz)','fXYHigh (Hz)', ...
%         'fSumHigh (Hz)','fLow (Hz)','Averaging windows'};
%     defaults = {'790','790','50000','25000','12000','100','200'};
%     param = getnumbers('Enter parameters:',prompts,defaults);
% end;

% k = 1.3806505e-23;
n = 9e-10; % [pN s/nm^2] water viscosity at 24.4C
% n = 9e-10*1.4; % [pN s/nm^2] water viscosity at 24.4C, for 10% glycerol
p = 1e-21; % [pN s^2/nm^4] water density
pbead = 1.05e-21; % [pN s^2/nm^4] bead density

beadA = param(1); gA = 3*pi*n*beadA; mA = pi/6*pbead*beadA^3;
% fsamp = param(3);
fXYHigh = param(4);
fSumHigh = param(5);
fLow = param(6);
avWin = param(7);

% Set Hydrodynamic Parameters
fmA = gA/(2*pi*(mA+pi/12*p*beadA^3)); fvA = n/p/(pi*beadA^2/4);

%can determine from timing data
fsamp = 1/data.sampperiod;

%hardwire background subtract off, doesn't work as is
% choice = 0;
% if choice == 1
%     data2 = ReadHighFreqFile(fsamp);
%     cd ..
%     cd 'Matlab Power Spectra\'
% 
%     AXSpecBkgd = AverageSpectrum(data2.A_X - mean(data2.A_X), fsamp, avWin);
%     AYSpecBkgd = AverageSpectrum(data2.A_Y - mean(data2.A_Y), fsamp, avWin);
%     ASumSpecBkgd = AverageSpectrum(data2.A_Sum - mean(data2.A_Sum), fsamp, avWin);
%     
%     cd ..
%     cd 'Analysis'
% else
    AXSpecBkgd = 0;
    AYSpecBkgd = 0;
    ASumSpecBkgd = 0;
% end;

% Change to Power Spectra Directory
olddir = pwd;
cd(analysisPath);

fit.path = startpath;
fit.file = file;

% Calculate Spectra
% display(['Calculating Power Spectra for ' file '...']);
[fit.AXSpec, fit.AXSpecRaw, fit.f] = AverageSpectrum(data.X - mean(data.X), fsamp, avWin);
[fit.AYSpec, fit.AYSpecRaw, fit.f] = AverageSpectrum(data.Y - mean(data.Y), fsamp, avWin);
[fit.ASumSpec, fit.ASumSpecRaw, fit.f] = AverageSpectrum(data.Sum - mean(data.Sum), fsamp, avWin);

% Subtract background
fit.AXSpec = fit.AXSpec - AXSpecBkgd;
fit.AYSpec = fit.AYSpec - AYSpecBkgd;
fit.ASumSpec = fit.ASumSpec - ASumSpecBkgd;

indXY = find(fit.f>fLow & fit.f<fXYHigh);
indSum = find(fit.f>fLow & fit.f<fSumHigh);
% display(['Done']);

% Use Analytical Lorentizian to Guess Initial Fit Values
[iguess_AXfc, iguess_AXD] = AnalyticalLorentzian(fit.AXSpec(indXY), fit.f(2), fLow, 10);
[iguess_AYfc, iguess_AYD] = AnalyticalLorentzian(fit.AYSpec(indXY), fit.f(2), fLow, 10);
[iguess_ASumfc, iguess_ASumD] = AnalyticalLorentzian(fit.ASumSpec(indSum), fit.f(2), fLow, 10);

% Fit with more sophisticated model
% display(['Fitting Spectra for ' file '...']);
switch lower(fitmodel)
    case {'hydrodynamicspectrum','hydro'}
        specFunc = @(f, para)HydrodynamicSpectrum(f, [para fvA fmA]);
        [fit.AXfc, fit.AXD] = NonLinearFit(fit.f(indXY), fit.AXSpec(indXY), specFunc, [iguess_AXfc iguess_AXD]);
        [fit.AYfc, fit.AYD] = NonLinearFit(fit.f(indXY), fit.AYSpec(indXY), specFunc, [iguess_AYfc iguess_AYD]);
        [fit.ASumfc, fit.ASumD] = NonLinearFit(fit.f(indSum), fit.ASumSpec(indSum), specFunc, [iguess_ASumfc iguess_ASumD]);
        predictedAX = HydrodynamicSpectrum(fit.f, [fit.AXfc fit.AXD fvA fmA]);
        predictedAY = HydrodynamicSpectrum(fit.f, [fit.AYfc fit.AYD fvA fmA]);
        predictedASum = HydrodynamicSpectrum(fit.f, [fit.ASumfc fit.ASumD fvA fmA]);

    case {'aliasedspectrum','aliased'}
        specFunc = @(f, para)AliasedLorentzianSpectrum(f, para, fsamp, 10);
        [fit.AXfc, fit.AXD] = NonLinearFit(fit.f(indXY), fit.AXSpec(indXY), specFunc, [iguess_AXfc iguess_AXD]);
        [fit.AYfc, fit.AYD] = NonLinearFit(fit.f(indXY), fit.AYSpec(indXY), specFunc, [iguess_AYfc iguess_AYD]);
        [fit.ASumfc, fit.ASumD] = NonLinearFit(fit.f(indSum), fit.ASumSpec(indSum), specFunc, [iguess_ASumfc iguess_ASumD]);
        predictedAX = AliasedLorentzianSpectrum(fit.f, [fit.AXfc fit.AXD], fsamp, 10);
        predictedAY = AliasedLorentzianSpectrum(fit.f, [fit.AYfc fit.AYD], fsamp, 10);
        predictedASum = AliasedLorentzianSpectrum(fit.f, [fit.ASumfc fit.ASumD], fsamp, 10);

    case {'aliasedhydrodynamicspectrum','aliasedhydro'}
        specFunc = @(f, para)AliasedHydrodynamicSpectrum(f, [para fvA fmA], fsamp, 10);
        [fit.AXfc, fit.AXD] = NonLinearFit(fit.f(indXY), fit.AXSpec(indXY), specFunc, [iguess_AXfc iguess_AXD]);
        [fit.AYfc, fit.AYD] = NonLinearFit(fit.f(indXY), fit.AYSpec(indXY), specFunc, [iguess_AYfc iguess_AYD]);
%         [fit.ASumfc, fit.ASumD] = NonLinearFit(fit.f(indSum), fit.ASumSpec(indSum), specFunc, [iguess_ASumfc iguess_ASumD]);
        predictedAX = AliasedHydrodynamicSpectrum(fit.f, [fit.AXfc fit.AXD fvA fmA],fsamp, 10);
        predictedAY = AliasedHydrodynamicSpectrum(fit.f, [fit.AYfc fit.AYD fvA fmA], fsamp, 10);
%         predictedASum = AliasedHydrodynamicSpectrum(fit.f, [fit.ASumfc fit.ASumD fvA fmA], fsamp, 10);

    case {'aliasedfilteredhydrodynamicspectrum','aliasedfiltered'} % To account for filtering of PSDs (don't use for QPD detectors!)
        calpara = CalParameters(mean(data.Sum), mean(data.Sum)) ;
        specFunc = @(f, para)AliasedFilteredHydrodynamicSpectrum(f, [para fvA fmA], fsamp, 10, [calpara.AX 0]);%[13e3 3.6 0.48 0]);%
        [fit.AXfc, fit.AXD] = NonLinearFit(fit.f(indXY), fit.AXSpec(indXY), specFunc, [iguess_AXfc iguess_AXD]);
        specFunc = @(f, para)AliasedFilteredHydrodynamicSpectrum(f, [para fvA fmA], fsamp, 10, [calpara.AX 0]);%[13e3 4.1 0.43 0]);%[calpara.AY 0]);%
        [fit.AYfc, fit.AYD] = NonLinearFit(fit.f(indXY), fit.AYSpec(indXY), specFunc, [iguess_AYfc iguess_AYD]);
        [fit.ASumfc, fit.ASumD] = NonLinearFit(fit.f(indSum), fit.ASumSpec(indSum), specFunc, [iguess_ASumfc iguess_ASumD]);
        predictedAX = AliasedFilteredHydrodynamicSpectrum(fit.f, [fit.AXfc fit.AXD fvA fmA],fsamp, 10,[calpara.AX 0]);%[13e3 3.6 0.48 0]);%[calpara.AX 0]); %
        predictedAY = AliasedFilteredHydrodynamicSpectrum(fit.f, [fit.AYfc fit.AYD fvA fmA], fsamp, 10,[calpara.AX 0]);%[13e3 4.1 0.43 0]);%[calpara.AY 0]); %
        predictedASum = AliasedHydrodynamicSpectrum(fit.f, [fit.ASumfc fit.ASumD fvA fmA], fsamp, 10);

    otherwise
        fit.AXfc = initguess_AXfc;
        fit.AXD = initguess_AXD;
        fit.AYfc = initguess_AYfc;
        fit.AYD = initguess_AYD;
        fit.ASumfc = initguess_ASumfc;
        fit.ASumD = initguess_ASumD;

        predictedAX = fit.AXD/pi^2*1./(fit.AXfc^2+fit.f.^2);
        predictedAY = fit.AYD/pi^2*1./(fit.AYfc^2+fit.f.^2);
        predictedASum = fit.ASumD/pi^2*1./(fit.ASumfc^2+fit.f.^2);

%         display(['Done']);

end
fit.predictedAX = predictedAX;
fit.predictedAY = predictedAY;

calibration.path = startpath;
calibration.file = file;
calibration.stamp = now;
calibration.date = date;
calibration.beadA = beadA;
calibration.fsamp = fsamp;
calibration.fXYHigh = fXYHigh;
calibration.fSumHigh = fSumHigh;
calibration.fLow = fLow;

calibration.alphaAX = sqrt(4.14/(gA*fit.AXD));
calibration.alphaAY = sqrt(4.14/(gA*fit.AYD));
% calibration.alphaASum = sqrt(4.14/(gA*fit.ASumD));

calibration.kappaAX = 2*pi*gA*fit.AXfc;
calibration.kappaAY = 2*pi*gA*fit.AYfc;
% calibration.kappaASum = 2*pi*gA*fit.ASumfc;

calibration.AXoffset = mean(data.X);
calibration.AYoffset = mean(data.Y);
calibration.AZoffset = mean(data.Sum);

fResHigh = 1000;
indRes = find(fit.f>fLow & fit.f<fResHigh);
calibration.AXmeanRes = mean(fit.AXSpec(indRes)./predictedAX(indRes)-1);
calibration.AYmeanRes = mean(fit.AYSpec(indRes)./predictedAY(indRes)-1);
% calibration.ASummeanRes = mean(fit.ASumSpec(indRes)./predictedASum(indRes)-1);

calibration.AXstdRes = std(fit.AXSpec(indRes)./predictedAX(indRes)-1);
calibration.AYstdRes = std(fit.AYSpec(indRes)./predictedAY(indRes)-1);
% calibration.ASumstdRes = std(fit.ASumSpec(indRes)./predictedASum(indRes)-1);

cd(olddir)


