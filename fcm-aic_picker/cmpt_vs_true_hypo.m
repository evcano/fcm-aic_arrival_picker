%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     cmpt_vs_true_hypo.m: compares computed and true hypocenters. 

close all;
clear;
clc;

% PARAMETERS --------------------------------------------------------------
project_name = 'synthetic_3';
arr_snr_dir = '../synthetic_data/synthetics_3_arr_snr';

ref_arr_dir = '../synthetic_data/synthetics_arr_times';

%% DO NOT EDIT BELOW THIS LINE --------------------------------------------

% load project information
project_file = strcat(project_name,'.project.mat');
load(project_file);

% load hypocenter location results
hypo_results_file = strcat('../hypocenter_location/hypo_',project_name,'.mat');
load(hypo_results_file)

tdom = par.tdom;
n_events = length(events_info.dir);
pick_res_rms = NaN(n_events,1);
pick_snr_rms = NaN(n_events,1);

for id = 1:n_events
    % load computed picks
    file = [num2str(events_info.id{id}) '_' events_info.name{id} '_results.mat'];
    load(fullfile(par.event_results_dir{id},file))
    c_ppick = event_results.p_pick .* par.dt;
    c_spick = event_results.s_pick .* par.dt;
    clear event_results;
    
    % load reference picks
    file = [events_info.name{id} '_arrivaltimes.mat'];
    ref = load(fullfile(ref_arr_dir,file));
    r_ppick = ref.p_pick(:) * par.dt;
    r_spick = ref.s_pick(:) * par.dt;
    clear ref;
    
    % load arrivals snr
    file = [events_info.name{id} '_arr_snr.mat'];
    load(fullfile(arr_snr_dir,file))
    r_psnr = 20*log10(p_snr);
    r_ssnr = 20*log10(s_snr);
    clear p_snr; clear s_snr;
    
    pick_res_rms(id) = rms(([r_ppick;r_spick] - [c_ppick;c_spick])*1000,'omitnan');
    pick_snr_rms(id) = rms([r_psnr;r_ssnr]);
end

% compute true and computed hypocenter difference
true_hypo = hypo_results.true_hyp;
true_az = hypo_results.true_az;
cmpt_hypo = hypo_results.cmpt_hyp;
cmpt_az = hypo_results.cmpt_az;
no_picks = hypo_results.npicks;

missing_hyp = no_picks == 0;
res_hypo = true_hypo - cmpt_hypo;
res_baz = abs(true_az - cmpt_az);
res_rms = rms(res_hypo,2);

% print results
nres = res_hypo(:,1);
eres = res_hypo(:,2);
zres = res_hypo(:,3);

nanmean(nres(~missing_hyp))
nanstd(nres(~missing_hyp))

nanmean(eres(~missing_hyp))
nanstd(eres(~missing_hyp))

nanmean(zres(~missing_hyp))
nanstd(zres(~missing_hyp))

%% FIGURES  ---------------------------------------------------------------
fz = 10; % font size

% fig 1: pick rms vs location rms
figure(1);
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 20, 10, 10]);
scatter(pick_res_rms(~missing_hyp),res_rms(~missing_hyp),40,no_picks(~missing_hyp),'filled');
grid on; box on;
xlabel('Arrival picking RMSE','fontsize',fz);
ylabel('Location RMSE','fontsize',fz);
h = colorbar;
ylabel(h,'No. of receivers');


% fig 2: baz rms vs location rms
figure;
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10, 20, 10, 10]);
scatter(res_baz(~missing_hyp),res_rms(~missing_hyp),40,no_picks(~missing_hyp),'filled');
grid on; box on;
xlabel('Back-azimuth RMSE','fontsize',fz);
ylabel('Location RMSE','fontsize',fz);
h = colorbar;
ylabel(h,'No. of receivers');


% fig 3: east residual vs north residual
figure;
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [0, 10, 10, 10]);
scatter(eres(~missing_hyp),nres(~missing_hyp),40,'k','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3);
xlabel('Easting residual (m)','fontsize',fz);
ylabel('Northing residual (m)','fontsize',fz);
box on; grid on;


% fig 4: north residual vs depth residual
figure;
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [10, 10, 10, 10]);
scatter(nres(~missing_hyp),zres(~missing_hyp),40,'k','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3);
xlabel('Northing residual (m)','fontsize',fz);
ylabel('Depth residual (m)','fontsize',fz);
box on; grid on;


% fig 5: east residual vs depth residual
figure;
set(gcf, 'color', 'w', 'units', 'centimeters', 'OuterPosition', [20, 10, 10, 10]);
scatter(eres(~missing_hyp),zres(~missing_hyp),40,'k','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3);
xlabel('Easting residual (m)','fontsize',fz);
ylabel('Depth residual (m)','fontsize',fz);
box on; grid on;
