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

function process_project(project_name,option)
fprintf('=================================\n');
fprintf('   FCM-AIC WAVE ARRIVAL PICKER   \n');
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
elseif strcmp(option,'single phase')
    events_in_queue = events_info.id(strcmp(events_info.status(:),'U'));
    if isempty(events_in_queue)
        fprintf('- There are no unkown-arrival events.\n');
    else
        fprintf(['- Identifying ' num2str(length(events_in_queue))  ' single phase events.\n']);
        identify_single_phase_event(project_name,events_info, par);
        fprintf('- Done.\n');
    end
    return
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
    
    % pick arrival times (flags = 0 unknown, 1=relabeled as p, 2=relabeled as s, 3=p relabeled as s
    event_results = struct('p_window',NaN(n_receivers,2),'p_pick',NaN(n_receivers,1),...
        's_window',NaN(n_receivers,2),'s_pick',NaN(n_receivers,1),...
        'p_moveout',NaN(n_receivers,1),'p_moveout_ni',NaN,'p_moveout_coef',NaN(1,3),...
        's_moveout',NaN(n_receivers,1),'s_moveout_ni',NaN,'s_moveout_coef',NaN(1,3),...
        'u_window',NaN(n_receivers,2),'u_pick',NaN(n_receivers,1),'pick_flag',NaN(n_receivers,1),...
        'baz',NaN(n_receivers,1)); 
    for r = 1:n_receivers
        fprintf(['    - Receiver ' num2str(r) ' ...\n']);
        tr_no = [1 2 3] + (r - 1)*3;
        tic;
        event_results = wave_arrival_picker(event_waveforms.amp(tr_no,:),event_results,r,par);
        etime = toc;
        fprintf(['\t* Done in ' num2str(etime) ' seconds.\n']);
    end
    
    % save figure of preliminary results
    if par.save_fig == 1
        fprintf('    - Saving figure of preliminary results ...\n');
        f = figure('visible','off');
        
        plot_waveforms(event_waveforms.amp,par.dt,'p_pick',event_results.p_pick,...
            's_pick',event_results.s_pick,'u_pick',event_results.u_pick,...
            'p_window',event_results.p_window,'s_window',event_results.s_window,...
            'u_window',event_results.u_window);
        
        fig_file = [num2str(event_id) '_' events_info.name{event_id} '_preliminary.jpg'];
        print(fullfile(par.event_results_dir{event_id},fig_file),'-djpeg');
        close(f);
        fprintf('\t* Done.\n');
    end
    
    % arrival times validation
    fprintf('    - Validating arrivals ...\n')
    tic;
    event_results = verify_phase_labeling(event_results, par);
    etime = toc;
    fprintf(['\t* Done in ' num2str(etime) ' seconds.\n']);
    
    % save figure of final results
    if par.save_fig == 1
        fprintf('    - Saving figure of results ...\n');
        f = figure('visible','off');
        
        plot_waveforms(event_waveforms.amp,par.dt,'p_pick',event_results.p_pick,...
            's_pick',event_results.s_pick,'u_pick',event_results.u_pick,...
            'p_window',event_results.p_window,'s_window',event_results.s_window,...
            'u_window',event_results.u_window,...
            'p_moveout',event_results.p_moveout,'s_moveout',event_results.s_moveout);
        
        fig_file = [num2str(event_id) '_' events_info.name{event_id} '.jpg'];
        print(fullfile(par.event_results_dir{event_id},fig_file),'-djpeg');
        close(f);
        fprintf('\t* Done.\n');
    end
    
    % save results
    output_file = [num2str(event_id) '_' events_info.name{event_id} '_results.mat'];
    save(fullfile(par.event_results_dir{event_id},output_file),'event_results');
    events_info = change_event_status(project_name,events_info,event_id,'D');
    % label event as single-phase
    if sum(event_results.pick_flag == 0) == n_receivers
        events_info = change_event_status(project_name,events_info,event_id,'U');
    end
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