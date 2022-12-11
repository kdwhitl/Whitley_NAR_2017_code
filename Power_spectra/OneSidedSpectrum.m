%  OneSidedSpectrum(timeSeries, fsamp, window) produces a properly
%  normalized power spectrum , P(f), of timeSeries, sampled at fsamp, using the
%  window function given by the user.  The resulting spectrum is plotted on
%  a loglog axis.  

%  Jeffrey Moffitt
%  November 15, 2004
function [Spectrum, f] = OneSidedSpectrum(timeSeries, fsamp, window)

timeSeries = timeSeries - mean(timeSeries);

% New Code 6/29/06  Make timeSeries an even number of points
timeSeries = timeSeries(1:(floor(length(timeSeries)/2)*2));

% New Code 8/3/06
N = length(timeSeries);

w = window(N);

windowedSeries = w'.*timeSeries;
windowSum = N*sum(w.*w);
transform = fft(windowedSeries);

power = transform.*conj(transform);

oneSidedPower = zeros(1, N/2+1);

% New Code 8/3/06
oneSidedPower(1) = power(1);  % Zero frequency point
oneSidedPower(2:(N/2)) = 2*power(2:(N/2));  % Positive and negative frequencies
%oneSidedPower(2:(N/2)) = power(2:(N/2)) + power((N/2 + 2):N));
oneSidedPower(N/2+1) = 2*power(N/2+1);  % fNyquist

df = fsamp/N;
f = 0:df:fsamp/2; %corrected 8/3/06 YC. 

Spectrum = oneSidedPower/(df*windowSum);  % Integrating the power spectrum will return
                                          % the observed RMS.



