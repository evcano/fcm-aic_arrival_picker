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

function [features] = compute_trace_features(tr, dt, mean_window, ...
    ppsd_window, sta_window, lta_window, features_taper)
% mean of the absolute value of the amplitude
tr_mean = compute_tr_mean(abs(tr), mean_window);
tr_mean(features_taper) = [];
tr_mean2 = minmaxn(tr_mean);
tr_mean2 = tr_mean2 - mean(tr_mean2);

% peak power spectral density
tr_ppsd = compute_tr_ppsd(tr, ppsd_window, 1/dt);
tr_ppsd(features_taper) = [];
tr_ppsd2 = minmaxn(tr_ppsd);
tr_ppsd2 = tr_ppsd2 - mean(tr_ppsd2);

% short-term average over long-term average
tr_stalta = compute_tr_stalta(abs(tr), sta_window, lta_window);
tr_stalta(features_taper) = [];
% tr_stalta(tr_stalta>10) = 10; % cap ampltidues larger than 10
tr_stalta2 = minmaxn(tr_stalta);
tr_stalta2 = tr_stalta2 - mean(tr_stalta2);

features = [tr_mean2' tr_ppsd2' tr_stalta2'];
end


function [tr_mean] = compute_tr_mean(tr, window)
if iscolumn(tr) == 1    
    tr = tr';
end

n = length(tr(1,:));
half_window = window * 0.5;
tr_mean = zeros(1,n);

for i = (1 + half_window):(n - half_window)
    window_start = i - half_window;
    window_end = i + half_window;
    tr_mean(i) = sum(tr(1,window_start:window_end)) * 1/n;
end
end


function [tr_ppsd] = compute_tr_ppsd(tr, window, df)
[~, ~, ~, psd] = spectrogram(tr, window, window - 1, [], df);
tr_ppsd = max(psd);
tr_ppsd = [tr_ppsd zeros(1,window - 1)];
end


function [tr_stalta] = compute_tr_stalta(tr, sta_window, lta_window)
if iscolumn(tr) == 1
    tr = tr';
end

n = length(tr(1,:));
tr_stalta = zeros(1,n);

for i = (1 + lta_window):n-sta_window
    lta = mean(tr(1,(i - lta_window):i));
    sta = mean(tr(1,i:i+sta_window));
    tr_stalta(i) = sta / lta;
end
tr_stalta(isinf(tr_stalta)) = 0;
tr_stalta(isnan(tr_stalta)) = 0;
end