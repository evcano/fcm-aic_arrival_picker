%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     locate_hypocenters: locate hypocenters of the synthetic datasets using
%                         picks obtained from the fcm-aic picker.                       

close all;
clear;
clc;
addpath 1d_ray_tracing_codes

% PARAMETERS --------------------------------------------------------------

project_name = 'synthetic_3';
ref_arrivals_dir = '../synthetic_data/synthetics_arr_times';


%% DO NOT EDIT BELOW THIS LINE ---------------------------------------------

project_dir = '../fcm-aic_picker/';

load sources;
load stations;
load model;
load(strcat(project_dir,project_name,'.project.mat'));

vmod(:,1) = z(2:end);
vmod(:,2) = vp;
vmod(:,3) = vs;
rec(:,3) = -rec(:,3);
n_events = size(events_info,1);

hypo_results = struct('true_hyp',NaN(n_events,4),'ini_hyp',NaN(n_events,4),...
    'cmpt_hyp',NaN(n_events,4),'true_az',NaN(n_events,1),...
    'cmpt_az',NaN(n_events,1),'npicks',NaN(n_events,1));

% inversion parameters
max_iter = 100;
stopping_criterion = 1e-9;
dx = 1; % jacobian inversion grid
dz = 1;

for id = 1:n_events
    % true hypocenter
    eno = str2double(events_info.name{id}(7:end));
    true_x = src_xyz(eno,:);
    true_x(3) = -true_x(3);
    
    % true azimuth
    azm = atan2d((true_x(2)-rec(1,2)),(true_x(1)-rec(1,1)));
    azm(azm<0) = 360 + azm(azm<0);
    true_azm = azm;
    
    % initial hypocenter guess (use true hypocenter and add some noise)
    x0 = [true_x(1); true_x(2); true_x(3); 0];
    noise = x0 .* [0.2; 0.2; 0.2; 0] .* (2*randi([0 1],[4,1])-1);
    x0 = x0 + noise;
    if x0(3) > vmod(end,1)
        x0(3) = vmod(end,1)-10;
    end
    x0(x0<0) = 0;
    
    % load picks
    file = [num2str(events_info.id{id}) '_' events_info.name{id} '_results.mat'];
    load(fullfile([project_dir par.event_results_dir{id}(2:end)],file))
    p_pick(:,1) = event_results.p_pick * par.dt;
    s_pick(:,1) = event_results.s_pick * par.dt;
    
    % estimated azimuth
    azm = nanmean(event_results.baz);
    azm(azm<0) = 360 + azm(azm<0);
    cmpt_azm = azm;
    
    % missing picks
    p_nan = ~isnan(p_pick);
    s_nan = ~isnan(s_pick);
    
    % observed data vector
    obT = [p_pick(p_nan); s_pick(s_nan)];
    
    % weight matrix
    W = eye(length(obT)+sum(p_nan));
    W(1:length(obT),:) = W(1:length(obT),:) * (1/0.002);
    W(length(obT)+1:end,:) = W(length(obT)+1:end,:) * (1/20);
    
    % observed data vector
    obT = [obT; azm.*ones(sum(p_nan),1)];
    
    % hypo inversion
    [cmpt_x,res,flag] = damped_ls_inv_hypoloc_3d(vmod, rec, obT, p_nan, s_nan, W, x0, dx, dz, max_iter, stopping_criterion);
    
    % comparison with true hypocenter
    clc;
    fprintf('\n True hypocenter: X=%0.3f, Y=%0.3f, Z=%0.3f \n',true_x(1),true_x(2),true_x(3))
    fprintf('\n Initial hypocenter res: X=%0.3f, Y=%0.3f, Z=%0.3f \n',true_x(1)-x0(1),true_x(2)-x0(2),true_x(3)-x0(3))
    fprintf('\n Comp hypocenter res: X=%0.3f, Y=%0.3f, Z=%0.3f \n',true_x(1)-cmpt_x(1),true_x(2)-cmpt_x(2),true_x(3)-cmpt_x(3))
    fprintf('\n Azm res = %0.3f \n',true_azm-cmpt_azm);
    
    hypo_results.true_hyp(id,:) = [true_x(:);0];
    hypo_results.true_az(id) = true_azm;
    hypo_results.ini_hyp(id,:) = x0(:);
    hypo_results.cmpt_hyp(id,:) = cmpt_x(:);
    hypo_results.cmpt_az(id) = cmpt_azm;
    hypo_results.npicks(id) = sum(p_nan);
end

filename = strcat('hypo_',project_name,'.mat');
save(filename, 'hypo_results');
