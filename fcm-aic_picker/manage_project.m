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

function manage_project(project_name, option)
fprintf('=======================\n');
fprintf('    PROJECT MANAGER    \n');
fprintf('=======================\n');
addpath('./src');

if strcmp(option,'create')
    fprintf(['\nCreating project "' project_name '":\n']);
    
    fprintf('\t- Reading parameters file ...\n');
    par = read_par_file(project_name);
    
    fprintf('\t- Validating parameters ...\n');
    validate_parameters(par);
    
    fprintf(['\t- Scanning events in "' par.datadir '" ...\n']);
    events_info = scan_events(par.datadir);
    
    fprintf('\t- Making neccesary directories ...\n');
    par = make_directories(project_name,par,events_info);
    
    fprintf('\t- Saving project file ...\n\n');
    save([project_name '.project.mat'],'par','events_info');
    
    fprintf(['Project "' project_name '" created succesfully.\n\n']);
elseif strcmp(option,'update')
    project_file = [project_name '.project.mat'];
    if exist(project_file,'file')
        fprintf(['\nUpdating project "' project_name '" parameters:\n']);
        load(project_file,'par');
        tmp = par.event_results_dir;
        
        fprintf('\t- Reading parameters file ...\n');
        par = read_par_file(project_name);
        
        fprintf('\t- Validating parameters ...\n');
        validate_parameters(par);
        par.event_results_dir = tmp;
        
        fprintf('\t- Saving project file ...\n\n');
        save([project_name '.project.mat'],'par','-append');
        
        fprintf(['Parameters of project "' project_name '" updated succesfully.\n\n']);
    else
        msg = ['The project "' project_name '" does not exist.'];
        error(msg);
    end
else
    msg = ['Unknow option "' option '". Use "create" to create a new project or "update"'...
        ' to update the parameters of an existing project.'];
    error(msg);
end
end


function [par] = read_par_file(project_name)
par_file = [project_name '.' 'par'];
if ~exist(par_file,'file')
    msg = ['Parameters file "' par_file '" does not exist.'];
    error(msg);
else
    file_content = fileread(par_file);
    file_parameters = regexp(file_content, '(?<field>\w+)=(?<value>[^\n]+)', 'names');
    par = struct();
    
    for i = 1:length(file_parameters)
        parameter = file_parameters(i).field;
        value = file_parameters(i).value;
        tmp = ismember(value,char([34 39]));
        value = strtrim(value(~tmp));
        value_tmp = str2double(value);
        if ~isnan(value_tmp)
            value = value_tmp;
        end
        par.(parameter) = value;
    end
end

% other parameters; do not change
par.n_clusters = 2;
par.fuzzifier = 2;
par.rectilinearity_threshold = 0.1;
end


function validate_parameters(par)
parameters_list = {'datadir';'data_format';'dt';'tdom';'mean_window';'ppsd_window';'sta_window';...
    'lta_window';'n_iterations';'stop_criteria';'save_fig'};

for i = 1:length(parameters_list)
    if ~isfield(par,parameters_list{i})
        msg = ['Missing parameter "' parameters_list{i} '".'];
        error(msg);
    end
end

if ~ischar(par.datadir) || exist(par.datadir,'dir') ~= 7
    msg = ['The data directory "' par.datadir '" doest not exist.'];
    error(msg);
end

if ~ischar(par.data_format) || ~strcmp(par.data_format,'SAC') && ~strcmp(par.data_format,'MAT') ...
        && ~strcmp(par.data_format,'sac') && ~strcmp(par.data_format,'mat')
    msg = ['Incorrect data format "' par.data_format '".'];
    error(msg);
end

if ischar(par.dt) || par.dt <= 0
    msg = 'Incorrect parameter "dt". It must be a positive number.';
    error(msg);
end

if ischar(par.tdom) || floor(par.tdom) ~= par.tdom || par.tdom <= 0
    msg = 'Incorrect parameter "tdom". It must be a positive integer.';
    error(msg);
end

if ischar(par.mean_window) || mod(par.mean_window,2) ~= 0 || par.mean_window < 2
    msg = 'Incorrect parameter "mean_window". It must be a positive even number.';
    error(msg);
end

if ischar(par.ppsd_window) || mod(par.ppsd_window,2) ~= 0 || par.ppsd_window < 2
    msg = 'Incorrect parameter "ppsd_window". It must be a positive even number.';
    error(msg);
end

if ischar(par.sta_window) || mod(par.sta_window,2) ~= 0 || par.sta_window < 2
    msg = 'Incorrect parameter "sta_window". It must be a positive even number.';
    error(msg);
end

if ischar(par.lta_window) || mod(par.lta_window,2) ~= 0 || par.lta_window < 2
    msg = 'Incorrect parameter "lta_window". It must be a positive even number.';
    error(msg);
end

if par.sta_window >= par.lta_window
    msg = 'Incorrect STA/LTA setting. "sta_window" must be smaller than "lta_window".';
    error(msg);
end

if ischar(par.n_iterations) || floor(par.n_iterations) ~= par.n_iterations || par.n_iterations < 1
    msg = 'Incorrect parameter "n_iterations". It must be a positive integer.';
    error(msg);
end

if ischar(par.stop_criteria) || par.stop_criteria <= 0
    msg = 'Incorrect parameter "stop_criteria". It must be a positive number.';
    error(msg);
end

if ischar(par.save_fig) || ~ismember(par.save_fig,[0 1 2])
    msg = 'Incorrect parameter "save_fig". Accepted values are "0","1", and "2".';
    error(msg);
end
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
    tmp.id{i,1} = i;
    tmp.name{i,1} = events_list{i};
    tmp.dir{i,1} = [datadir '/' events_list{i}];
    tmp.status{i,1} = 'P';
end

events_info = table(tmp.id,tmp.name,tmp.dir,tmp.status,'VariableName',{'id','name','dir','status'});
end

function par = make_directories(project_name,par,events_info)
results_dir = ['./' project_name '_results'];
if exist(results_dir,'dir') == 0
    mkdir(results_dir);
end

for i = 1:height(events_info)
    event_results_dir = [results_dir '/' num2str(events_info.id{i}) '_' events_info.name{i}];
    par.event_results_dir{i} = event_results_dir;
    if exist(event_results_dir,'dir') == 0
        mkdir(event_results_dir);
    end
end
end