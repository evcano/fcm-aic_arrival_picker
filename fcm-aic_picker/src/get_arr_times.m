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

function [p_pick, s_pick] = get_arr_times(waveforms,p_window,s_window,tdom,win_ext_factor,plot_results)
% extend p window
p_window(1) = p_window(1) - win_ext_factor*tdom;

% compute p wave SNR
waveforms_noise = waveforms(:,1:(p_window(1) - tdom));
p_wave = waveforms(:,p_window(1):p_window(2));
p_wave_snr = rms(p_wave,2) ./ rms(waveforms_noise,2);

% identify the component with highest p wave SNR
p_wave_snr(2) = 0;
[~,p_component] = max(p_wave_snr);

% pick p wave arrival
p_pick = pick_arrival_1(waveforms(p_component,:),p_window,plot_results);

if isempty(s_window)
    s_pick = NaN;
else
    % extend s wave window
    s_window(1) = max(p_pick + tdom, s_window(1) - win_ext_factor*tdom);
    
    % s wave will be picked on the components different from p_wave_component
    s_components = 1:3;
    s_components(s_components == p_component) = [];
    
    % pick s wave arrival
    s_pick = pick_arrival_2(waveforms(s_components,:),s_window,plot_results);
end
end


function aic = akaike_ic(waveform)
n_samples = length(waveform);
aic = zeros(1,n_samples-1);

for i = 1:(n_samples - 1)
    a = var(waveform(1:i));
    b = var(waveform(i + 1:n_samples));
    if a > 0
        a = log(a);
    else
        a = 0;
    end
    if b > 0
        b = log(b);
    else
        b = 0;
    end
    aic(i) = i * a + (n_samples - i - 1) * b;
end
end


function [pick,aic] = pick_arrival_1(waveform,window,plot_results)
arrival = waveform(window(1):window(2));
aic = akaike_ic(arrival);
aic = minmaxn(aic);
[~,pick_tmp] = min(aic);
pick = pick_tmp + window(1) - 1;

if plot_results
    figure(8);
    fig_size = [10, 10, 10, 10];
    set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
    subplot(211);
    plot(arrival);
    hold on;
    plot(pick_tmp,arrival(pick_tmp),'or');
    subplot(212);
    plot(aic);
    hold on;
    plot(pick_tmp,aic(pick_tmp),'or');
end
end


function [pick,aic] = pick_arrival_2(waveforms,window,plot_results)
n_waveforms = size(waveforms,1);
aic = zeros(n_waveforms,(window(2) - window(1)));
arrival = zeros(n_waveforms,(1 + window(2) - window(1)));

for i = 1:n_waveforms
    arrival(i,:) = waveforms(i,window(1):window(2));
    aic(i,:) = akaike_ic(arrival(i,:));
    aic(i,:) = minmaxn(aic(i,:));
end

aic_stack = sum(aic) .* (1 / n_waveforms);
[~,pick_tmp] = min(aic_stack);
pick = pick_tmp + window(1) - 1;

if plot_results
    figure(7)
    fig_size = [0, 12, 10, 15];
    set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
    subplot(511); plot(arrival(1,:)); hold on; plot(pick_tmp,arrival(1,pick_tmp),'or');
    subplot(512); plot(aic(1,:));
    subplot(513); plot(arrival(2,:)); hold on; plot(pick_tmp,arrival(2,pick_tmp),'or');
    subplot(514); plot(aic(2,:));
    subplot(515); plot(aic_stack); hold on; plot(pick_tmp,aic_stack(pick_tmp),'or');
end
end
