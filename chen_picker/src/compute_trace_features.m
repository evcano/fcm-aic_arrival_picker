%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


function [features] = compute_trace_features(tr, mean_window, ...
    pow_window, sta_window, lta_window, features_taper)
% mean
tr_mean = compute_tr_mean(abs(tr), mean_window);
tr_mean(features_taper) = [];
tr_mean2 = minmaxn(tr_mean);
tr_mean2 = tr_mean2 - mean(tr_mean2);

% power
tr_pow = compute_tr_pow(abs(tr), pow_window);
tr_pow(features_taper) = [];
tr_pow2 = minmaxn(tr_pow);
tr_pow2 = tr_pow2 - mean(tr_pow2);

% short-term average over long-term average
tr_stalta = compute_tr_stalta(abs(tr), sta_window, lta_window);
tr_stalta(features_taper) = [];
tr_stalta2 = minmaxn(tr_stalta);
tr_stalta2 = tr_stalta2 - mean(tr_stalta2);

features = [tr_mean2' tr_pow2' tr_stalta2'];
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


function [tr_pow] = compute_tr_pow(tr, window)
if iscolumn(tr) == 1    
    tr = tr';
end

n = length(tr(1,:));
half_window = window * 0.5;
tr_pow = zeros(1,n);

for i = (1 + half_window):(n - half_window)
    window_start = i - half_window;
    window_end = i + half_window;
    tr_pow(i) = sum( tr(1,window_start:window_end).^2 );
end
end


function [tr_stalta] = compute_tr_stalta(tr, sta_window, lta_window)
if iscolumn(tr) == 1
    tr = tr';
end

n = length(tr(1,:));
tr_stalta = zeros(1,n);

for i = (1 + lta_window):n-sta_window
    lta = mean(tr(1,(i-lta_window):i));
    sta = mean(tr(1,i-sta_window:i));
    tr_stalta(i) = sta / lta;
end
tr_stalta(isinf(tr_stalta)) = 0;
tr_stalta(isnan(tr_stalta)) = 0;
end