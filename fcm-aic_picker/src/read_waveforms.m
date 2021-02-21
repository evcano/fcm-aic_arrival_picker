%     FCM-AIC WAVE ARRIVAL PICKER
%     ---------------------------
%     Copyright (C) November 2020  Eduardo Valero Cano,
%     King Abdullah University of Science and Technology (KAUST).
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function waveforms = read_waveforms(input, file_format)
input_type = exist(input);

if input_type == 7
    waveforms_list = dir([input '/*.' file_format]);
    waveforms_list = {waveforms_list.name}';
    n_waveforms = length(waveforms_list);
    if n_waveforms == 0
        msg = ['There are no "' file_format '" files in "' input '".'];
        error(msg);
    else
        for i = 1:n_waveforms
            waveforms_list{i} = [input '/' waveforms_list{i}];
        end
    end
elseif input_type == 2
    if endsWith(input,file_format)
        waveforms_list{1} = input;
        n_waveforms = 1;
    else
        msg = ['The format of "' input '" is not "' file_format '".'];
        error(msg);
    end
else
    msg = ['The directory/file "' input '" does not exist.'];
    error(msg);
end

if strcmp(file_format,'SAC') || strcmp(file_format,'sac')
    for i = 1:n_waveforms
        tmp = rsac(waveforms_list{i});
        waveforms.time(i,:) = tmp(:,1);
        waveforms.amp(i,:) = tmp(:,2);
        waveforms.header(i,:) = tmp(:,3);
        waveforms.file(i,:) = waveforms_list{i};
    end
elseif strcmp(file_format,'MAT') || strcmp(file_format,'mat')
    for i = 1:n_waveforms
        tmp = load(waveforms_list{i});
        waveforms.time(i,:) = tmp.waveform(:,1);
        waveforms.amp(i,:) = tmp.waveform(:,2);
        waveforms.header(i,:) = tmp.waveform(:,3);
        waveforms.file(i,:) = waveforms_list{i};
    end
else
    msg = ['Incorrect file format "' file_format '". Accepted formats: "SAC","MAT".'];
    error(msg);
end
end


function [varargout] = rsac(varargin)
%RSAC    Read SAC binary files.
%    RSAC('sacfile') reads in a SAC (seismic analysis code) binary
%    format file into a 3-column vector.
%    Column 1 contains time values.
%    Column 2 contains amplitude values.
%    Column 3 contains all SAC header information.
%    Default byte order is big-endian.  M-file can be set to default
%    little-endian byte order.
%
%    usage:  output = rsac('sacfile')
%
%    Examples:
%
%    KATH = rsac('KATH.R');
%    plot(KATH(:,1),KATH(:,2))
%
%    [SQRL, AAK] = rsac('SQRL.R','AAK.R');
%
%    by Michael Thorne (4/2004)   mthorne@asu.edu

for nrecs = 1:nargin
    
    sacfile = varargin{nrecs};
    
    %---------------------------------------------------------------------------
    %    Default byte-order
    %    endian  = 'big-endian' byte order (e.g., UNIX)
    %            = 'little-endian' byte order (e.g., LINUX)
    
    endian = 'little-endian';
    
    if strcmp(endian,'big-endian')
        fid = fopen(sacfile,'r','ieee-be');
    elseif strcmp(endian,'little-endian')
        fid = fopen(sacfile,'r','ieee-le');
    end
    
    % read in single precision real header variables:
    %---------------------------------------------------------------------------
    for i=1:70
        h(i) = fread(fid,1,'single');
    end
    
    % read in single precision integer header variables:
    %---------------------------------------------------------------------------
    for i=71:105
        h(i) = fread(fid,1,'int32');
    end
    
    
    % Check header version = 6 and issue warning
    %---------------------------------------------------------------------------
    % If the header version is not NVHDR == 6 then the sacfile is likely of the
    % opposite byte order.  This will give h(77) some ridiculously large
    % number.  NVHDR can also be 4 or 5.  In this case it is an old SAC file
    % and rsac cannot read this file in.  To correct, read the SAC file into
    % the newest verson of SAC and w over.
    %
    if (h(77) == 4 || h(77) == 5)
        message = strcat('NVHDR = 4 or 5. File: "',sacfile,'" may be from an old version of SAC.');
        error(message)
    elseif h(77) ~= 6
        message = strcat('Current rsac byte order: "',endian,'". File: "',sacfile,'" may be of opposite byte-order.');
        error(message)
    end
    
    % read in logical header variables
    %---------------------------------------------------------------------------
    for i=106:110
        h(i) = fread(fid,1,'int32');
    end
    
    % read in character header variables
    %---------------------------------------------------------------------------
    for i=111:302
        h(i) = (fread(fid,1,'char'))';
    end
    
    % read in amplitudes
    %---------------------------------------------------------------------------
    
    YARRAY     = fread(fid,'single');
    
    if h(106) == 1
        XARRAY = (linspace(h(6),h(7),h(80)))';
    else
        error('LEVEN must = 1; SAC file not evenly spaced')
    end
    
    % add header signature for testing files for SAC format
    %---------------------------------------------------------------------------
    h(303) = 77;
    h(304) = 73;
    h(305) = 75;
    h(306) = 69;
    
    % arrange output files
    %---------------------------------------------------------------------------
    OUTPUT(:,1) = XARRAY;
    OUTPUT(:,2) = YARRAY;
    OUTPUT(1:306,3) = h(1:306)';
    
    %pad xarray and yarray with NaN if smaller than header field
    if h(80) < 306
        OUTPUT((h(80)+1):306,1) = NaN;
        OUTPUT((h(80)+1):306,2) = NaN;
    end
    
    fclose(fid);
    
    varargout{nrecs} = OUTPUT;
end
end