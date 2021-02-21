%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     cmpt_vs_ref_arr_times.m: compares computed and reference arrival picks. 

close all;
clear;
clc;

% PARAMETERS --------------------------------------------------------------
project_name = 'synthetic_1';
arr_snr_dir = '../synthetic_data/synthetics_3_arr_snr';

ref_arr_dir = '../synthetic_data/synthetics_arr_times';

%% DO NOT EDIT BELOW THIS LINE --------------------------------------------

% load project information
project_file = strcat(project_name,'.project.mat');
load(project_file);

tdom = par.tdom;
n_events = length(events_info.dir);

% initialize variables
c_ppick = [];
c_spick = [];
c_upick = [];
c_pick_flag = [];
r_ppick = [];
r_spick = [];
pwave_snr = [];
swave_snr = [];

for id = 1:n_events
    % load computed picks
    file = [num2str(events_info.id{id}) '_' events_info.name{id} '_results.mat'];
    load(fullfile(par.event_results_dir{id},file))
    c_ppick = [c_ppick; event_results.p_pick .* par.dt];
    c_spick = [c_spick; event_results.s_pick .* par.dt];
    c_upick = [c_upick; event_results.u_pick .* par.dt];
    c_pick_flag = [c_pick_flag; event_results.pick_flag];
    clear event_results;
    
    % load reference picks
    file = [events_info.name{id} '_arrivaltimes.mat'];
    ref = load(fullfile(ref_arr_dir,file));
    r_ppick = [r_ppick; ref.p_pick(:) * par.dt];
    r_spick = [r_spick; ref.s_pick(:) * par.dt];
    clear ref;
    
    % load arrivals snr
    file = [events_info.name{id} '_arr_snr.mat'];
    load(fullfile(arr_snr_dir,file))
    pwave_snr = [pwave_snr; p_snr];
    swave_snr = [swave_snr; s_snr];
    clear p_snr; clear s_snr;
end

% convert arrivals snr to dB
pwave_snr = (20*log10(pwave_snr));
swave_snr = (20*log10(swave_snr));

% pick residuals (in miliseconds)
p_res = (r_ppick - c_ppick) * 1000;
s_res = (r_spick - c_spick) * 1000;

% missing computed picks
mis_c_ppick = isnan(c_ppick);
mis_c_spick = isnan(c_spick);

% relabeled picks
rel_c_ppick = c_pick_flag == 1;
rel_c_spick = c_pick_flag == 2 | c_pick_flag==3;

% missing reference picks
mis_r_ppick = isnan(r_ppick);
mis_r_spick = isnan(r_spick);

% keep residuals between + - threshold
thr = 50;
tmp = find(p_res >= -thr & p_res <= thr);
p_res2 = p_res(tmp);

% keep residuals between + - threshold
thr = 50;
tmp = find(s_res >= -thr & s_res <= thr);
s_res2 = s_res(tmp);

% print statistics
fprintf('P picks residuals mean: %5.3f ms\n',nanmean(p_res2));
fprintf('P picks residuals std: %5.3f ms\n',nanstd(p_res2));
fprintf('No. of P picks residuals: %i \n',numel(p_res2));
fprintf('\nS picks residuals mean: %5.3f ms\n',nanmean(s_res2));
fprintf('S picks residuals std: %5.3f ms\n',nanstd(s_res2));
fprintf('No. of S picks residuals: %i \n',numel(s_res2));

% save missed picks
filename = strcat('mis_c_ppick_',project_name,'.mat');
save(filename, 'mis_c_ppick')

%% FIGURES  ---------------------------------------------------------------
fz = 10; % figures font size

% fig 1: missing, picked, and relabeled picks p picks
figure(1)
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 20, 10, 10])
X = [sum(mis_c_ppick), sum(~isnan(c_ppick))-sum(rel_c_ppick), sum(rel_c_ppick)];
p = pie(X,[1 0 0]);
p(1).FaceColor = [0.6 0.3 0.0];
p(2).String = [p(2).String ' (' num2str(X(1)) ')'];
p(3).FaceColor = [0 0.5 1];
p(4).String = [p(4).String ' (' num2str(X(2)) ')'];
p(5).FaceColor = 'g';
p(6).String = [p(6).String ' (' num2str(X(3)) ')'];
legend({'Missing','Picked','Relabeled'})


% fig 2: missing, picked, and relabeled picks s picks
figure(2)
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10 20, 10, 10])
X = [sum(mis_c_spick), sum(~isnan(c_spick))-sum(rel_c_spick), sum(rel_c_spick)];
p = pie(X,[1 0 0]);
p(1).FaceColor = [0.6 0.3 0.0];
p(2).String = [p(2).String ' (' num2str(X(1)) ')'];
p(3).FaceColor = [0 0.5 1];
p(4).String = [p(4).String ' (' num2str(X(2)) ')'];
p(5).FaceColor = 'g';
p(6).String = [p(6).String ' (' num2str(X(3)) ')'];
legend({'Missing','Picked','Relabeled'})


% fig 3: p pick residual vs snr
exclude = mis_c_ppick + rel_c_ppick;
figure(3);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 10, 10, 10])
scatter(p_res(~exclude),pwave_snr(~exclude),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.3);
hold on;
scatter(p_res(rel_c_ppick),pwave_snr(rel_c_ppick),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.3);
plot([10 10], [min(pwave_snr) max(pwave_snr)], 'k--', 'linewidth', 1.5)
plot([-10 -10], [min(pwave_snr) max(pwave_snr)], 'k--', 'linewidth', 1.5)
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Signal-to-noise ratio (dB)', 'fontsize', fz)
grid on; 
box on; 
axis('tight');


% fig 4: s pick residual vs snr
exclude = mis_c_spick + rel_c_spick;
figure(4);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10, 10, 10, 10])
scatter(s_res(~exclude),swave_snr(~exclude),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r','MarkerFaceAlpha',0.3);
hold on;
scatter(s_res(rel_c_spick),swave_snr(rel_c_spick),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.3);
plot([10 10], [min(swave_snr) max(swave_snr)], 'k--', 'linewidth', 1.5)
plot([-10 -10], [min(swave_snr) max(swave_snr)], 'k--', 'linewidth', 1.5)
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Signal-to-noise ratio (dB)', 'fontsize', fz)
grid on; 
box on; 
axis('tight');

% fig 5: p pick residual histogram
edges = min(p_res):5:max(p_res);
figure(5);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 0, 10, 10])
sh = histogram(p_res,'FaceColor',[0, 0, 1],'FaceAlpha',1,'EdgeColor',[0 0 0],'EdgeAlpha',1,'BinEdges',edges);
sd = fitdist(p_res2,'normal');
ys = pdf(sd,edges);
ys = ys./max(ys)*max(sh.Values);
hold on;
plot(edges, ys,'k','linewidth',1);
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Counts', 'fontsize', fz)
grid on; box on;

% fig 6: s pick residual histogram
edges = min(s_res):5:max(s_res);
figure(6);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10, 0, 10, 10])
sh = histogram(s_res,'FaceColor',[1, 0, 0],'FaceAlpha',1,'EdgeColor',[0 0 0],'EdgeAlpha',1,'BinEdges',edges);
sd = fitdist(s_res2,'normal');
ys = pdf(sd,edges);
ys = ys./max(ys)*max(sh.Values);
hold on;
plot(edges, ys,'k','linewidth',1);
xlabel('Residual (ms)', 'fontsize', fz)
ylabel('Counts', 'fontsize', fz)
grid on; box on;
