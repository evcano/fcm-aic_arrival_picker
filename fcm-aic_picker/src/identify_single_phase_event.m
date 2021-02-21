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

function identify_single_phase_event(project_name,events_info, par)
n_events = size(events_info,1);

% RANSAC parameters
sampleSize = par.ransac_samplesize;
maxDistance = par.ransac_maxdistance;
fitFunction = @(data) polyfit(data(:,1),data(:,2),2);
misfitFunction = @(model, data) sum((data(:,2) - polyval(model,data(:,1))).^2,2);

% load p and s wave moveouts
for e = 1:n_events
    event_id = e;
    input_file = [num2str(event_id) '_' events_info.name{event_id} '_results.mat'];
    load(fullfile(par.event_results_dir{event_id},input_file));
    
    p_moveout(:,e) = (event_results.p_moveout - min(event_results.p_moveout)) * par.dt;
    p_moveout_ni(e) = event_results.p_moveout_ni;
    
    s_moveout(:,e) = (event_results.s_moveout - min(event_results.s_moveout)) * par.dt;
    s_moveout_ni(e) = event_results.s_moveout_ni;
end

events_in_queue = events_info.id(strcmp(events_info.status(:),'U'));
for e = 1:length(events_in_queue)
    event_id = events_in_queue{e};
    input_file = [num2str(event_id) '_' events_info.name{event_id} '_results.mat'];
    load(fullfile(par.event_results_dir{event_id},input_file));
    
    u_picks = event_results.u_pick;
    rno(:,1) = 1:length(u_picks);
    x(:,1) = u_picks;
    tmp = isnan(x);
    modelRANSAC = ransac([rno(~tmp),x(~tmp)], fitFunction, misfitFunction,...
        sampleSize, maxDistance);
    u_moveout = round(polyval(modelRANSAC,rno));
    u_moveout = (u_moveout - min(u_moveout)) * par.dt;
    dif = abs(u_moveout - u_picks*par.dt);
    inliers = find(dif <= maxDistance);
    outliers = find(dif > maxDistance);
    
    figure;
    fig_size = [10, 0, 7, 17];
    fz = 10;
    set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
    plot(p_moveout(:,p_moveout_ni>10),[1:20],'color',[0 0 1 0.5],'LineWidth',0.5);
    hold on;
    plot(s_moveout(:,s_moveout_ni>10),[1:20],'color',[1 0 0 0.5],'LineWidth',0.5);
    plot(u_moveout,[1:20],'g','LineWidth',2);
    box on; grid on; axis('tight');
    xlabel('Time (s)','fontsize',fz);
    ylabel('Receiver','fontsize',fz);
    set(gca,'ydir','reverse','ytick',1:20);
    close all;
    
    pdif = nanmean(rms(p_moveout - u_moveout,1).^2);
    sdif = nanmean(rms(s_moveout - u_moveout,1).^2);
    if sdif < pdif
        tmp1 = inliers;
        tmp2 = outliers;
    else
        tmp1 = outliers;
        tmp2 = inliers;
    end
    
    event_results.s_pick(tmp1) = event_results.u_pick(tmp1);
    event_results.s_window(tmp1,:) = event_results.u_window(tmp1,:);
    event_results.pick_flag(tmp1) = 2;
    
    event_results.p_pick(tmp2) = event_results.u_pick(tmp2);
    event_results.p_window(tmp2,:) = event_results.u_window(tmp2,:);
    event_results.pick_flag(tmp2) = 1;
    
    event_results.u_pick = NaN(size(event_results.u_pick));
    event_results.u_window = NaN(size(event_results.u_window));
    
    output_file = [num2str(event_id) '_' events_info.name{event_id} '_results.mat'];
    save(fullfile(par.event_results_dir{event_id},output_file),'event_results');
    events_info = change_event_status(project_name,events_info,event_id,'DU');
end
end

function events_info = change_event_status(project_name,events_info,id,status)
events_info.status{id} = status;
save([project_name '.project.mat'],'events_info','-append');
end