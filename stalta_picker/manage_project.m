%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

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
end


function validate_parameters(par)
parameters_list = {'datadir';'data_format';'dt';'tdom';'sta_window';...
    'lta_window';'stalta_thr';'save_fig'};

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