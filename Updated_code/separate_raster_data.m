function data = separate_raster_data(data,smooth_scan)
%just takes data file and divides data into raster parts

N = data.nscanfb * 2;%total number of scans (forward + backward)

if nargin == 1
    smooth_scan = 1;
end

%after 090415, sometimes may be extra data at the end of set after scans
%have finished.  Check for this and separate it.  Almost certain after
%090903.
if str2double(data.file(1:6)) > 90414
    %this assumes that the raw sampling rate is once per trap cycle,
    %but right now I'm running 1 out of 3 (slower rate for 66 kHz
    %timesharing, 100211).  But the data file doesn't have a way of
    %figuring this out.  So cheat.  Do something special if rawperiod =
    %6e-5.
    if data.rawperiod < 7e-5 && data.rawperiod > 5e-5
        NRasterpts = floor(N*data.nsteps*data.ncycles/(data.sampperiod/data.rawperiod)/4);
    else
        NRasterpts = floor(N*data.nsteps*data.ncycles/(data.sampperiod/data.rawperiod));
    end
    
    %090909, mjc, after changing labview code to give more general and
    %physical scan parameters (scan rate, size, step size are in nm no
    %counts or MHz), the save data integration rate doesn't always
    %exactly line up with scan length.  Try this fix (throws out a few
    %extra points).
    NRasterpts = floor(NRasterpts/N)*N;
    
    if smooth_scan == 0
        %this way will have equal number of points per step
        NRasterpts = floor(NRasterpts/(N*data.nsteps))*N*data.nsteps;
    end
    
    Ntpts = length(data.time);
    if Ntpts > NRasterpts
        %separate extra data from scans
        data.A_Xafter = data.A_X(NRasterpts*data.nchannels+1:end);
        data.A_Yafter = data.A_Y(NRasterpts*data.nchannels+1:end);
        data.B_Xafter = data.B_X(NRasterpts*data.nchannels+1:end);
        data.B_Yafter = data.B_Y(NRasterpts*data.nchannels+1:end);
        data.C_Xafter = data.C_X(NRasterpts*data.nchannels+1:end);
        data.C_Yafter = data.C_Y(NRasterpts*data.nchannels+1:end);
        data.A_Sumafter = data.A_Sum(NRasterpts*data.nchannels+1:end);
        data.B_Sumafter = data.B_Sum(NRasterpts*data.nchannels+1:end);
        
        data.A_X = data.A_X(1:NRasterpts);
        data.A_Y = data.A_Y(1:NRasterpts);
        data.B_X = data.B_X(1:NRasterpts);
        data.B_Y = data.B_Y(1:NRasterpts);
        data.C_X = data.C_X(1:NRasterpts);
        data.C_Y = data.C_Y(1:NRasterpts);
        data.A_Sum = data.A_Sum(1:NRasterpts);
        data.B_Sum = data.B_Sum(1:NRasterpts);
        
        if data.chanpattern == 2
            data.A_FB_Xafter = data.A_FB_X(NRasterpts+1:end);
            data.B_FB_Xafter = data.B_FB_X(NRasterpts+1:end);
            
            data.A_FB_X = data.A_FB_X(1:NRasterpts);
            data.B_FB_X = data.B_FB_X(1:NRasterpts);
        elseif data.chanpattern == 4
            data.A_FB_Yafter = data.A_FB_Y(NRasterpts+1:end);
            data.A_FB_Xafter = data.A_FB_X(NRasterpts+1:end);
            data.B_FB_Yafter = data.B_FB_Y(NRasterpts+1:end);
            data.B_FB_Xafter = data.B_FB_X(NRasterpts+1:end);
            
            data.A_FB_Y = data.A_FB_Y(1:NRasterpts);
            data.A_FB_X = data.A_FB_X(1:NRasterpts);
            data.B_FB_Y = data.B_FB_Y(1:NRasterpts);
            data.B_FB_X = data.B_FB_X(1:NRasterpts);
        end
    end
else %no extra data for sure
    NRastertpts = length(data.time);
end

data.NRasterpts = NRasterpts;

%separate data into scans
data.samples_scan = NRasterpts/N;%number of samples per one scan

% data.samples_scan
% length(data.time)

data.t = data.time(1:data.samples_scan);

%reshape arrays to separate out scans
data.A_X = reshape(data.A_X,[],N)';
data.A_Y = reshape(data.A_Y,[],N)';
data.B_X = reshape(data.B_X,[],N)';
data.B_Y = reshape(data.B_Y,[],N)';
data.C_X = reshape(data.C_X,[],N)';
data.C_Y = reshape(data.C_Y,[],N)';
data.A_Sum = reshape(data.A_Sum,[],N)';
data.B_Sum = reshape(data.B_Sum,[],N)';

if data.chanpattern == 2
    data.A_FB_X = reshape(data.A_FB_X,[],N)';
    data.B_FB_X = reshape(data.B_FB_X,[],N)';
elseif data.chanpattern == 4
    data.A_FB_Y = reshape(data.A_FB_Y,[],N)';
    data.A_FB_X = reshape(data.A_FB_X,[],N)';
    data.B_FB_Y = reshape(data.B_FB_Y,[],N)';
    data.B_FB_X = reshape(data.B_FB_X,[],N)';
end

%flip the even rows (the backwards scans)
b = (1:data.nscanfb)*2;
for j = b
    data.A_X(j,:) = fliplr(data.A_X(j,:));
    data.A_Y(j,:) = fliplr(data.A_Y(j,:));
    data.B_X(j,:) = fliplr(data.B_X(j,:));
    data.B_Y(j,:) = fliplr(data.B_Y(j,:));
    data.C_X(j,:) = fliplr(data.C_X(j,:));
    data.C_Y(j,:) = fliplr(data.C_Y(j,:));
    data.A_Sum(j,:) = fliplr(data.A_Sum(j,:));
    data.B_Sum(j,:) = fliplr(data.B_Sum(j,:));
    
    if data.chanpattern == 2
        data.A_FB_X(j,:) = fliplr(data.A_FB_X(j,:));
        data.B_FB_X(j,:) = fliplr(data.B_FB_X(j,:));
    elseif data.chanpattern == 4
        data.A_FB_Y(j,:) = fliplr(data.A_FB_Y(j,:));
        data.A_FB_X(j,:) = fliplr(data.A_FB_X(j,:));
        data.B_FB_Y(j,:) = fliplr(data.B_FB_Y(j,:));
        data.B_FB_X(j,:) = fliplr(data.B_FB_X(j,:));
    end
end

%moving trap position
%set of scan frequency values
if smooth_scan == 1 %compute assuming a smooth linear scan
    data.nu = (0:(data.samples_scan-1))/(data.samples_scan-1)*data.scanrange;%change in moving trap freq.
else %compute assuming actual steps
    data.nu_reduced = (0:(data.nsteps-1))/(data.nsteps-1)*data.scanrange;
    nsampstep = data.samples_scan/data.nsteps;%how many samples we have per step
    reptemp = repmat(data.nu_reduced,nsampstep,1);%rep the array by how many samples
    data.nu = reshape(reptemp,[],data.samples_scan);%now a row vector
end

data.nudelta = (data.t2x-data.t1x) - data.nu;

