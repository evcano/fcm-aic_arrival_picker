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

function event_results = wave_arrival_picker(waveforms,event_results,receiver_no,par)
% create taper for trace features
n_samples = length(waveforms);
left_taper = 1:par.lta_window;
right_taper = (n_samples - par.mean_window + 1):n_samples;
features_taper = [left_taper, right_taper];

features = cell(3,1);
signal_membership = zeros(3,n_samples-length(features_taper));

% clustering of receiver components using fuzzy c-means
for i = 1:3
    features{i} = compute_trace_features(waveforms(i,:),par.dt,par.mean_window,...
        par.ppsd_window,par.sta_window,par.lta_window,features_taper);
    
    memberships = fuzzy_c_means(features{i},par.n_clusters,par.n_iterations,par.fuzzifier,...
        par.stop_criteria);
    
    signal_membership(i,:) = determine_signal_cluster(memberships);
end

% signal membership stacking
s_signal_membership = sum(signal_membership,1) .* (1 / 3);

signal_threshold = par.signal_threshold * mean(s_signal_membership);
s_signal_membership = [zeros(1,length(left_taper)), s_signal_membership, zeros(1,length(right_taper))];

% candidate-arrivals windows computation
candidate_arr_windows = get_candidate_arr_windows(s_signal_membership,signal_threshold,par.tdom);

if isempty(candidate_arr_windows)
    fprintf('\t* No candidate arrival windows identified.\n');
    event_results = end_arrival_picker(event_results,receiver_no);
    return;
end

% arrival windows computation
[p_window,s_window,u_window,waveforms_rot,baz] = get_arr_windows(waveforms,candidate_arr_windows,par.rectilinearity_threshold);

% arrival time picking
if ~isnan(u_window) % unknow wave arrival
    fprintf('\t* Unknow arrival identified.\n');
    event_results.pick_flag(receiver_no) = 0;
    u_pick = get_arr_times(waveforms_rot,u_window,[],par.tdom,par.win_ext,par.plot_results);
    p_pick = NaN;
    s_pick = NaN;
else % p and s wave arrival
    fprintf('\t* P-wave and S-wave arrival identified.\n');
    u_pick = NaN;
    [p_pick,s_pick] = get_arr_times(waveforms_rot,p_window,s_window,par.tdom,par.win_ext,par.plot_results);
end

% output results
event_results.p_window(receiver_no,:) =  p_window;
event_results.p_pick(receiver_no,:) =  p_pick;
event_results.s_window(receiver_no,:) =  s_window;
event_results.s_pick(receiver_no,:) =  s_pick;
event_results.u_window(receiver_no,:) =  u_window;
event_results.u_pick(receiver_no,:) =  u_pick;
event_results.baz(receiver_no,:) = baz;

if par.plot_results
    plot_results(waveforms,receiver_no, features, features_taper,...
        signal_membership, s_signal_membership, signal_threshold, event_results, candidate_arr_windows)
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


function event_results = end_arrival_picker(event_results,receiver_no)
fprintf('\t* Arrival picking failed.\n');
event_results.p_window(receiver_no,:) =  NaN(1,2);
event_results.p_pick(receiver_no,:) =  NaN;
event_results.s_window(receiver_no,:) =  NaN(1,2);
event_results.s_pick(receiver_no,:) =  NaN;
event_results.u_window(receiver_no,:) =  NaN(1,2);
event_results.u_pick(receiver_no,:) =  NaN;
end


function plot_results(waveforms,receiver_no, features, features_taper, signal_membership, s_signal_membership,...
    signal_threshold, event_results,candidate_arr_windows)
t = 1:length(waveforms);
tt = t; tt(features_taper) = [];
waveformst = waveforms; waveformst(:,features_taper) = [];

figure(9);
fig_size = [0, 0, 10, 10];
set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
for ii = 1:3
    subplot(5,3,ii);
    plot(tt,waveformst(ii,:),'k');
    subplot(5,3,ii+3);
    plot(tt,features{ii}(:,1),'k');
    subplot(5,3,ii+6);
    plot(tt,features{ii}(:,2),'k');
    subplot(5,3,ii+9);
    plot(tt,features{ii}(:,3),'k');
    subplot(5,3,ii+12);
    plot(tt,signal_membership(ii,:),'k');
end

%
figure(10);
fig_size = [10, 0, 10, 10];
set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
for ii = 1:3
    subplot(4,1,ii);
    plot_window(waveforms(ii,:),event_results.p_window(receiver_no,:),[0 0 1 0.2],10);
    hold on;
    plot_window(waveforms(ii,:),event_results.s_window(receiver_no,:),[1 0 0 0.2],10);
    plot_window(waveforms(ii,:),event_results.u_window(receiver_no,:),[0 1 0 0.2],10);
    for jj = 1:size(candidate_arr_windows,1)
        if sum(candidate_arr_windows(jj,:)-event_results.p_window(receiver_no,:)) ~=0 && ...
                sum(candidate_arr_windows(jj,:)-event_results.s_window(receiver_no,:)) ~=0 ...
                && sum(candidate_arr_windows(jj,:)-event_results.u_window(receiver_no,:)) ~=0
            
            plot_window(waveforms(ii,:),candidate_arr_windows(jj,:),[0 0 0 0.2],10);
        end
    end
    plot(t,waveforms(ii,:),'k');
    plot_pick(t,waveforms(ii,:),event_results.p_pick(receiver_no),'b','--',10);
    plot_pick(t,waveforms(ii,:),event_results.s_pick(receiver_no),'r','--',10);
    axis('tight')
    box on;
end

%
subplot(4,1,4);
plot_window(s_signal_membership,event_results.p_window(receiver_no,:),[0 0 1 0.2],10);
plot_window(s_signal_membership,event_results.s_window(receiver_no,:),[1 0 0 0.2],10);
plot_window(s_signal_membership,event_results.u_window(receiver_no,:),[0 1 0 0.2],10);
for jj = 1:size(candidate_arr_windows,1)
    if sum(candidate_arr_windows(jj,:)-event_results.p_window(receiver_no,:)) ~=0 && ...
            sum(candidate_arr_windows(jj,:)-event_results.s_window(receiver_no,:)) ~=0 ...
            && sum(candidate_arr_windows(jj,:)-event_results.u_window(receiver_no,:)) ~=0
        
        plot_window(s_signal_membership,candidate_arr_windows(jj,:),[0 0 0 0.2],10);
    end
end
hold on;
plot(t,s_signal_membership,'k');
plot([t(1) t(end)],[signal_threshold signal_threshold],'r');
axis('tight')
box on;

fz=10;

fig_size = [10, 0, 7, 4];
figure(20)
set(gcf,'color','w','units','centimeters','OuterPosition',fig_size);
plot_window(s_signal_membership,candidate_arr_windows(jj,:)*0.0005,[0 1 0 0.2],20);
hold on;
plot(t*0.0005,s_signal_membership,'k');
box on; grid on; axis('tight')
set(gca,'fontsize',fz);
xlabel('Time (s)','fontsize',fz);

pause;
close all;
end


function plot_pick(time,waveform,pick,color,marker,figno)
if ~isnan(pick)
    figure(figno)
    plot([time(pick) time(pick)],[min(waveform) max(waveform)],'Color',color,'LineStyle',marker);
end
end


function plot_window(waveform,window,color,figno)
if ~isnan(window)
    figure(figno)
    rectangle('Position',[window(1) min(waveform) (window(2)-window(1)) (max(waveform)-min(waveform))],...
        'FaceColor',color,'EdgeColor',color);
end
end