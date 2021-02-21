%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

function event_results = wave_arrival_picker(waveforms,event_results,tr_no,par)
sta_window = par.sta_window;
lta_window = par.lta_window;
thr = par.stalta_thr;
left_taper = 1:par.lta_window;

for i = 1:3
    stalta = compute_tr_stalta(abs(waveforms(i,:)), sta_window, lta_window);
    stalta(left_taper) = [];
    stalta = stalta / max(stalta);
    [~,tmp] = findpeaks(stalta,'MinPeakHeight',thr); % arrival = (first sta/lta peak > thr)
    
    event_results.p_pick(tr_no(i)) = tmp(1) + par.lta_window;
end
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