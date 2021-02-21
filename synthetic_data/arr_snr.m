%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     arr_snr.m: computes p- and s-wave SNR using the reference picks.

close all;
clear;
clc;


% PARAMETERS --------------------------------------------------------------
ref_arrivals_dir = './synthetics_arr_times';
data_dir = './synthetics_3';

data_format = 'MAT';
out_dir = './synthetics_arr_snr';
tdom = 50; % in samples


% DO NOT EDIT BELOW THIS LINE ---------------------------------------------
events_info = scan_events(data_dir);
n_events = length(events_info.dir);

for e = 1:n_events
    % load event
    event_waveforms = read_waveforms(events_info.dir{e}, data_format);
    n_receivers = size(event_waveforms.amp, 1) / 3;
    
    % load reference picks
    file = [events_info.name{e} '_arrivaltimes.mat'];
    tmp = load(fullfile(ref_arrivals_dir, file));
    ref_p_pick = tmp.p_pick;
    ref_s_pick = tmp.s_pick;
    
    % compute snr
    p_snr = nan(n_receivers,1);
    s_snr = nan(n_receivers,1);
    
    for r = 1:n_receivers
        c  = [1 2 3] + 3*(r-1);
        
        % p arrival snr
        if isnan(ref_p_pick(r)) == 0
            noise_win(1) = 1;
            noise_win(2) = ref_p_pick(r) - tdom;
            
            noise = event_waveforms.amp(c, noise_win(1):noise_win(2));
            p_wave = event_waveforms.amp(c, ref_p_pick(r):ref_p_pick(r) + tdom);
            
            noise_rms = rms(noise(:));
            p_wave_rms = rms(p_wave(:));
            p_snr(r) = mean(p_wave_rms ./ noise_rms);
        else
            p_snr(r) = NaN;
        end
        
        % s arrival snr
        if isnan(ref_s_pick(r)) == 0
            if isnan(ref_p_pick(r)) == 0
                noise_win(1) = 1;
                noise_win(2) = ref_p_pick(r) - tdom;
            else
                noise_win(1) = 1;
                noise_win(2) = ref_s_pick(r) - tdom;
            end
            
            noise = event_waveforms.amp(c, noise_win(1):noise_win(2));
            s_wave = event_waveforms.amp(c, ref_s_pick(r):ref_s_pick(r) + tdom);
            
            noise_rms = rms(noise(:));
            s_wave_rms = rms(s_wave(:));
            s_snr(r) = mean(s_wave_rms ./ noise_rms);
        else
            s_snr(r) = NaN;
        end
    end
    
    file = [events_info.name{e} '_arr_snr.mat'];
    save(fullfile(out_dir, file), 'p_snr', 's_snr')
end


function [events_info] = scan_events(datadir)
events_list = dir(datadir);
events_list = events_list([events_list(:).isdir]==1); % list only subdirectories
events_list = events_list(~ismember({events_list.name}, {'.','..'})); % delete "." and ".."
events_list = {events_list.name}';

n_events = length(events_list);
if n_events == 0
    msg = ['There are no events in "' datadir '".'];
    error(msg);
end

tmp = struct();
for i = 1:n_events
    tmp.name{i,1} = events_list{i};
    tmp.dir{i,1} = [datadir '/' events_list{i}];
end

events_info = table(tmp.name,tmp.dir,'VariableName',{'name','dir'});
end