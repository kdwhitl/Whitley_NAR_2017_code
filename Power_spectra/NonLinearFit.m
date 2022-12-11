function [fc, D, sfc, sD] = NonLinearFit(f, spec, specFunc, init);
% [fc, D, sfc, sD, mse] = NonLinearFit(f, spec, specFunc, init) returns the
% corner frequency, fc, and diffusion constant, D, from a nonlinear least
% squares fit to the function Pexp/Ptheory - 1.  This is the desired least
% squares function for spectra that have been averaged.  
% Ptheory is defined by specFunc and must be defined as specFunc(f, [fc,
% D]).  init are the inital guesses for fc and D, [fc, D].  

% Jeff Moffitt
% 8/5/06

% Define Pexp/Ptheory for the fitting routine
func = @(para, f)spec./specFunc(f, para);

% Call the nonlinear least squares minimization routine
[paraFit, resid, J] = nlinfit(f, ones(1, length(spec)), func, init);

% Assign fit values
fc = paraFit(1);
D = paraFit(2);

% Calculate 95% Confidence Intervals for the fit
ci = nlparci(real(paraFit), real(resid), real(J));  % Kludge!!

% Calculate 1 sigma from these intervals
sfc = (ci(1,2) - ci(1,1))/4;
sD = (ci(2,2) - ci(2,1))/4;

% Assign the mean squared error (or Chi square)
% mse = sigma(1,1);
% Not avaliable in all versions of Matlab ???