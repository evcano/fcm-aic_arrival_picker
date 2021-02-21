%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


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
    
    % plot waveforms
    plot(x, ye, 'k', 'LineWidth', 1);
    plot(x, yn, 'k', 'LineWidth', 1);
    plot(x, yz, 'k', 'LineWidth', 1);
    
    plot_pick(p_pick(t-2),ye,dt,'ob',1);
    plot_pick(p_pick(t-1),yn,dt,'ob',1);
    plot_pick(p_pick(t),yz,dt,'ob',1);
    
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


function plot_pick(pick,wave,dt,color,line_width)
if isnan(pick)
    return;
else
    plot(pick*dt,wave(pick),color,'LineWidth',line_width);
end
end