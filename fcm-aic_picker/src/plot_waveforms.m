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

function plot_waveforms(waveforms,varargin)
[n_waveforms, n_samples] = size(waveforms);
if n_waveforms > n_samples
    waveforms = waveforms';
end
[n_waveforms, n_samples] = size(waveforms);
n_receivers = n_waveforms / 3;
nan_array = nan(n_receivers,1);

p = inputParser;
addOptional(p,'dt',1);
addParameter(p,'p_window',nan_array);
addParameter(p,'p_pick',nan_array);
addParameter(p,'s_window',nan_array);
addParameter(p,'s_pick',nan_array);
addParameter(p,'u_window',nan_array);
addParameter(p,'u_pick',nan_array);
addParameter(p,'p_moveout',nan_array);
addParameter(p,'s_moveout',nan_array);
parse(p,varargin{:});
dt = p.Results.dt;
p_window = p.Results.p_window;
p_pick = p.Results.p_pick;
s_window = p.Results.s_window;
s_pick = p.Results.s_pick;
u_window = p.Results.u_window;
u_pick = p.Results.u_pick;
p_moveout = p.Results.p_moveout;
s_moveout = p.Results.s_moveout;

dx = 0.5;
weigths = dx ./ max(abs(waveforms),[],2);
waveforms = waveforms .* weigths; % normalize waveforms

x = (1:n_samples) .* dt;
ticky = cell(n_receivers * 3, 1);
hold on
for r = 1:n_receivers
    t = r * 3;
    ye = -waveforms(t - 2,:) + 2 * (t - 2) * dx; % change sign due to axis reversal
    yn = -waveforms(t - 1,:) + 2 * (t - 1) * dx;
    yz = -waveforms(t,:) + 2 * t * dx;
    
    % plot arrival windows
    plot_window(p_window(r,:),r,dt,[0 0 1 0.15]);
    plot_window(s_window(r,:),r,dt,[1 0 0 0.15]);
    plot_window(u_window(r,:),r,dt,[0 1 0 0.15]);
    
    % plot arrival picks
    plot_pick(p_pick(r),r,dt,'b',1.2);
    plot_pick(s_pick(r),r,dt,'r',1.2);
    plot_pick(u_pick(r),r,dt,'g',1.2);
    
    % moveout
    plot_pick(p_moveout(r),r,dt,'b--',1.3);
    plot_pick(s_moveout(r),r,dt,'r--',1.3);
    
    % plot waveforms
    plot(x, ye, 'k', 'LineWidth', 1);
    plot(x, yn, 'k', 'LineWidth', 1);
    plot(x, yz, 'k', 'LineWidth', 1);
    
    ticky{t - 2} = [num2str(r) ' - E'];
    ticky{t - 1} = ['N'];
    ticky{t} = ['Z'];
end

% figure appearance
fig_size = [0, 0, 1, 1];
fig_font_size = 12;
lab_font_size = 12;

set(gcf,'color','w','units','Normalized','OuterPosition',fig_size);
set(gca,'fontsize',fig_font_size);
set(gca,'ydir','reverse','ytick',1:n_waveforms,'yticklabel',{ticky{:}});
ylabel('Receiver','fontsize',lab_font_size);
if dt == 1
    xlabel('Sample No.','fontsize',lab_font_size);
else
    xlabel('Time (s)','fontsize',lab_font_size);
end
axis([0 n_samples*dt 0 dx*(2*n_waveforms+1)]);
grid on;
box on;
hold off;
end


function plot_pick(pick,r,dt,color,line_width)
if isnan(pick)
    return;
else
    plot([pick*dt pick*dt],[r*3-2.5 r*3+0.5],color,'LineWidth',line_width);
end
end


function plot_window(window,r,dt,color)
if isnan(window)
    return;
else
    rectangle('Position',[window(1)*dt, r*3-2.5, (window(2)*dt - window(1)*dt) 3],...
        'FaceColor',color,'EdgeColor',color);
end
end