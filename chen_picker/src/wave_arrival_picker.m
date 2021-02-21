%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


function event_results = wave_arrival_picker(waveforms,event_results,tr_no,par)
% create taper for trace features
n_samples = length(waveforms);
left_taper = 1:par.lta_window;
if par.mean_window > par.pow_window
    right_taper = (n_samples - par.mean_window + 1):n_samples;
else
    right_taper = (n_samples - par.pow_window + 1):n_samples;
end
features_taper = [left_taper, right_taper];

features = cell(3,1);
signal_membership = zeros(3,n_samples-length(features_taper));

% clustering of receiver components using fuzzy c-means
for i = 1:3
    features{i} = compute_trace_features(waveforms(i,:),par.mean_window,...
        par.pow_window,par.sta_window,par.lta_window,features_taper);
    
    memberships = fuzzy_c_means(features{i},par.n_clusters,par.n_iterations,par.fuzzifier,...
        par.stop_criteria);
    
    signal_membership(i,:) = determine_signal_cluster(memberships);
    
    % taper signal membership
    signal_membership(i,1:2*par.tdom) = 0;
    signal_membership(i,end-2*par.tdom:end) = 0;
    
    thr = par.pick_thr * mean(signal_membership(i,:));
    tmp = find(signal_membership(i,:) > thr);
    
    event_results.p_pick(tr_no(i)) = tmp(1) + left_taper(end);
end
end


function [signal_membership] = determine_signal_cluster(memberships)
c1_mean_membership = mean(memberships(:,1));
c2_mean_membership = mean(memberships(:,2));

if c1_mean_membership < c2_mean_membership
    signal_membership = memberships(:,1);
else
    signal_membership = memberships(:,2);
end
end