%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


function process_project(project_name,option)
fprintf('=================================\n');
fprintf('   CHEN (2020) ARRIVAL PICKER   \n');
fprintf('=================================\n');
addpath('./src');

fprintf('\n- Reading project file ...\n');
[par,events_info] = read_project_file(project_name);

% get events in queue
if strcmp(option,'all')
    tmp_flag = 0;
    while tmp_flag == 0
        msg = ['- Do you want to process all the events in the project "' project_name ...
            '"? This will delete any previous results (yes/no):'];
        inp = input(msg,'s');
        
        if strcmp(inp,'yes') || strcmp(inp,'y')
            events_in_queue = events_info.id(:);
            tmp_flag = 1;
        elseif strcmp(inp,'no') || strcmp(inp,'n')
            fprintf('- Operation cancelled.\n');
            return
        else
            fprintf('- Incorrect answer. Enter "yes" or "no".\n');
            tmp_flag = 0;
        end
    end
elseif strcmp(option,'pending')
    events_in_queue = events_info.id(strcmp(events_info.status(:),'P'));
    if isempty(events_in_queue)
        fprintf(['- There are no pending events. Use "all" to process all events or '...
            'enter the ID of an event to process.\n']);
        return
    end
elseif isnumeric(option)
    if ismember(option,cell2mat(events_info.id))
        events_in_queue = {option};
    else
        msg = ['There is no event with ID "' num2str(option) '".'];
        error(msg);
    end
else
    msg = ['Incorrect option. Use "pending" to process pending events; "all" to process all events; or enter '...
        'the ID of an event to process.'];
    error(msg);
end

n_events = length(events_in_queue);
fprintf(['- A total of ' num2str(n_events) ' events will be processed.\n\n']);

% events in queue loop
for e = 1:n_events
    % load waveforms
    event_id = events_in_queue{e};
    event_dir = events_info.dir{event_id};
    fprintf(['- Picking arrivals of event "' events_info.name{event_id} '" (ID = ' num2str(event_id) ') ...\n']);
    events_info = change_event_status(project_name,events_info,event_id,'P');
    event_waveforms = read_waveforms(event_dir,par.data_format);
    n_receivers = size(event_waveforms.amp,1) / 3;
    
    % pick arrival times
    event_results = struct('p_pick',NaN(n_receivers*3,1));
    for r = 1:n_receivers
        fprintf(['    - Receiver ' num2str(r) ' ...\n']);
        tr_no = [1 2 3] + (r - 1)*3;
        tic;
        event_results = wave_arrival_picker(event_waveforms.amp(tr_no,:),event_results,tr_no,par);
        etime = toc;
        fprintf(['\t* Done in ' num2str(etime) ' seconds.\n']);
    end
    
    % save figure of final results
    if par.save_fig == 1
        fprintf('    - Saving figure of results ...\n');
        f = figure('visible','off');
        
        plot_waveforms(event_waveforms.amp,par.dt,'p_pick',event_results.p_pick);
        
        fig_file = [num2str(event_id) '_' events_info.name{event_id} '.jpg'];
        print(fullfile(par.event_results_dir{event_id},fig_file),'-djpeg');
        close(f);
        fprintf('\t* Done.\n');
    end
    
    % save results
    output_file = [num2str(event_id) '_' events_info.name{event_id} '_results.mat'];
    save(fullfile(par.event_results_dir{event_id},output_file),'event_results');
    events_info = change_event_status(project_name,events_info,event_id,'D');
end
end


function [par,events_info] = read_project_file(project_name)
project_file = [project_name '.project.mat'];
if ~exist(project_file,'file')
    msg = ['Project file "' project_file '"does not exist.'];
    error(msg);
else
    load(project_file,'par','events_info');
    fprintf(['- Project file "' project_file '" read succesfully.\n']);
end
end


function events_info = change_event_status(project_name,events_info,id,status)
events_info.status{id} = status;
save([project_name '.project.mat'],'events_info','-append');
end