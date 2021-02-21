%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     cmpt_vs_ref_arr_times.m: compares computed and reference arrival picks. 

clc;
clear;
close all;


% PARAMETERS --------------------------------------------------------------
project_name = 'synthetic_1';
arr_snr_dir = '../synthetic_data/synthetics_3_arr_snr';
mis_ppicks_file = '../fcm-aic_picker/mis_c_ppick_synthetic_3.mat';

ref_arrivals_dir = '../synthetic_data/synthetics_arr_times';

%% DO NOT EDIT BELOW THIS LINE --------------------------------------------

% load project information
project_file = strcat(project_name,'.project.mat');
load(project_file);

tdom = par.tdom;
n_events = length(events_info.dir);

c_ppick = [];
r_ppick = [];
pwave_snr = [];

for id = 1:n_events
    % load computed picks
    file = [num2str(events_info.id{id}) '_' events_info.name{id} '_results.mat'];
    load(fullfile(par.event_results_dir{id},file))
    
    % load reference picks
    file = [events_info.name{id} '_arrivaltimes.mat'];
    ref = load(fullfile(ref_arrivals_dir,file));
    r_ppick = [r_ppick; ref.p_pick(:) * par.dt];
    
    % load arrivals snr
    file = [events_info.name{id} '_arr_snr.mat'];
    load(fullfile(arr_snr_dir,file))
    pwave_snr = [pwave_snr; p_snr];
    clear p_snr;
    
    % determine pick with minimum residual
    for i = 1:20
        trno = [1 2 3] + 3*(i-1);
        tmp = abs(event_results.p_pick(trno) - ref.p_pick(i)).*par.dt;
        [~,ind] = min(tmp);
        c_ppick_tmp(i,1) = event_results.p_pick(trno(ind)).* par.dt;
    end
    c_ppick = [c_ppick; c_ppick_tmp];
    
    % save best pick per receiver
    filename = [num2str(id) '_' events_info.name{id} '_best_results.mat'];
    save(fullfile(par.event_results_dir{id},filename),'c_ppick_tmp');
    
    clear event_results;
end

% pick residual (in miliseconds)
p_res = (r_ppick - c_ppick) * 1000;

% keep residuals between + - threshold
thr = 50;
tmp = find(p_res >= -thr & p_res <= thr);
p_res2 = p_res(tmp);

% convert snr to dB
pwave_snr = (20*log10(pwave_snr));

% load p picks missed by fcm-aic
load(mis_ppicks_file)
exclude = mis_c_ppick;

% print statistics
fprintf('P picks residuals mean: %5.3f ms\n',nanmean(p_res2));
fprintf('P picks residuals std: %5.3f ms\n',nanstd(p_res2));
fprintf('No. of P picks residuals: %i \n',numel(p_res2));

%% FIGURES  ---------------------------------------------------------------
fz = 10; % font size

% fig 1: p pick residual vs snr
figure(1);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 20, 10, 10])
scatter(p_res(~exclude),pwave_snr(~exclude),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3);
hold on;
scatter(p_res(mis_c_ppick),pwave_snr(mis_c_ppick),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.6 0.3 0.0],'MarkerFaceAlpha',0.3);
plot([10 10], [min(pwave_snr) max(pwave_snr)], 'k--', 'linewidth', 1.5)
plot([-10 -10], [min(pwave_snr) max(pwave_snr)], 'k--', 'linewidth', 1.5)
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Signal-to-noise ratio (dB)', 'fontsize', fz)
grid on; 
box on; 
axis('tight');

% fig 2: p pick residual histogram
edges = min(p_res):5:max(p_res);
figure(2);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10, 20, 10, 10])
sh = histogram(p_res,'FaceColor',[0, 0, 1],'FaceAlpha',1,'EdgeColor',[0 0 0],'EdgeAlpha',1,'BinEdges',edges);
sd = fitdist(p_res2,'normal');
ys = pdf(sd,edges);
ys = ys./max(ys)*max(sh.Values);
hold on;
plot(edges, ys,'k','linewidth',1);
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Counts', 'fontsize', fz)
grid on; 
box on;
