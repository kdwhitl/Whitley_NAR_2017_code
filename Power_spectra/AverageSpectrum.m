%  AverageSpectrum(timeSeries, fsamp, numWindows, window) produces a properly
%  normalized power spectrum , P(f), of timeSeries, sampled at fsamp, using the
%  window function given by the user, and averaged over numWindows of windows.
%  The resulting spectrum is plotted on a loglog axis.

%  Jeffrey Moffitt
%  November 15, 2004
%
% Modified to allow smoothing of data: simple smoothing of outlier data
% points
% 081215, mjc
%
% Incorporated tiny nowindow function. 170309 kw

function [spectrum,spectrumraw, f] = AverageSpectrum(timeSeries, fsamp, numWindows, window)

% RemOutlier = 0;%turn remove outlier feature on/off, 1/0

%default parameter values
if nargin < 4
%     window = @nowindow;
    window = @(num) ones(num,1)/num;
end

if nargin < 3
    numWindows = 1;
%     window = @nowindow;
    window = @(num) ones(num,1)/num;
end

index1 = 1;
index2 = floor(length(timeSeries)/numWindows);
[sumSpectrum, f] = OneSidedSpectrum(timeSeries(index1:index2), fsamp, window);
[~, c] = size(sumSpectrum);
tempSpectrum = zeros(numWindows,c);
tempSpectrum(1,:) = sumSpectrum;
for i = 2:numWindows
    index1 = floor(length(timeSeries)/numWindows)*(i-1)+1;
    index2 = floor(length(timeSeries)/numWindows)*i;
    tempSpectrum(i,:) = OneSidedSpectrum(timeSeries(index1:index2), fsamp, window);
end
sumSpectrum = sum(tempSpectrum);
spectrum = sumSpectrum/numWindows;
%loglog(f, spectrum);

spectrumraw = spectrum;

% if RemOutlier
% Sigma = std(spectrum);
% Mean = mean(spectrum);

%     outliers = (spectrum - Mean) > 2*Sigma;
%     spectrum2 = spectrum;
%     
%     hfig = figure;
%     subplot(2,1,1)
%     loglog(f,spectrum);
%     title('raw data')
%     subplot(2,1,2)
%     temp = circshift(spectrum,[0,-5]);
%     spectrum2(outliers) = temp(outliers);
%     loglog(f,spectrum2);
%     title('2 sigma outliers removed')
%     spectrum = spectrum2;
% end


end


