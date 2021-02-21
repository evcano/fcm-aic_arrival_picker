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

function [p_window, s_window, u_window, tr_rotated, baz] = get_arr_windows(tr, candidate_arr_windows,...
    rectilinearity_threshold)
n_candidate_arr_windows = size(candidate_arr_windows, 1);

% first case: one candidate arrival = unknown arrival
if n_candidate_arr_windows == 1
    tr_rotated = tr;
    u_window(1,1:2) = candidate_arr_windows;
    p_window(1,1:2) = NaN;
    s_window(1,1:2) = NaN;
    baz = NaN;
elseif n_candidate_arr_windows > 1
    [E1,~,E3] = compute_candidate_arr_windows_eigenvalues(tr,candidate_arr_windows);
    candidate_arr_windows_rectilinearity =  1 - E3 ./ E1;
    [p_window, p_window_index] = get_p_window(candidate_arr_windows,candidate_arr_windows_rectilinearity,rectilinearity_threshold);
    
    if p_window_index == n_candidate_arr_windows
        % second case: candidate arrival max rectilinearity = unknown arrival
        tr_rotated = tr;
        u_window = p_window;
        p_window(1,1:2) = NaN;
        s_window(1,1:2) = NaN;
        baz = NaN;
    else
        % third case: p and s arrivals
        [tr_rotated, baz] = rotate_ray_centered(tr,p_window);
        u_window(1,1:2) =  NaN;
        s_window = get_s_window(tr_rotated, candidate_arr_windows, p_window_index);
    end
end
end


% function to compute eigen values
function [E1,E2,E3] = compute_candidate_arr_windows_eigenvalues(tr, candidate_arr_windows)
n_candidate_arr_windows = size(candidate_arr_windows, 1);

E1 = zeros(n_candidate_arr_windows,1);
E2 = zeros(n_candidate_arr_windows,1);
E3 = zeros(n_candidate_arr_windows,1);

for i = 1:n_candidate_arr_windows
    [~,S,~] = svd(tr(:,candidate_arr_windows(i,1):candidate_arr_windows(i,2))');
    tmp = diag(S).^2;
    E1(i) = tmp(1);
    E2(i) = tmp(2);
    E3(i) = tmp(3);
end
end


function [p_window, p_window_index] = get_p_window(candidate_arr_windows,candidate_arr_windows_rectilinearity,rectilinearity_threshold)
max_rectilinearity = max(candidate_arr_windows_rectilinearity);
candidate_arr_windows_rectilinearity((max_rectilinearity - candidate_arr_windows_rectilinearity) > rectilinearity_threshold*max_rectilinearity) = 0;

p_window_index = find(candidate_arr_windows_rectilinearity > 0);
p_window_index = p_window_index(1);
p_window = candidate_arr_windows(p_window_index(1),:);
end


function [tr_rotated, baz1] = rotate_ray_centered(tr, p_window)
% 180 degree search for backazimuth
p_wave_E = tr(1,p_window(1):p_window(2));
p_wave_N = tr(2,p_window(1):p_window(2));
p_wave_energy_EN = zeros(180,1);
for i = 1:180
    p_wave_energy_EN(i) = sum((-p_wave_E * cosd(i) + p_wave_N * sind(i)).^2);
end

% backazimuth minimizes p wave energy on the East and North components
[~, baz1] = min(p_wave_energy_EN);

tr_rtv(1,:) = (-tr(1,:) * sind(baz1)) - (tr(2,:) * cosd(baz1)); % radial
tr_rtv(2,:) = (-tr(1,:) * cosd(baz1)) + (tr(2,:) * sind(baz1)); % transverse
tr_rtv(3,:) = tr(3,:); % vertical

p_wave_R = tr_rtv(1,p_window(1):p_window(2));
p_wave_V = tr_rtv(3,p_window(1):p_window(2));
p_wave_energy_RV = zeros(90,1);
for i = 1:90
    p_wave_energy_RV(i) = sum((p_wave_R * cosd(i) + p_wave_V * sind(i)).^2);
end

[~, baz2] = max(p_wave_energy_RV);

R = [cosd(baz2) 0 sind(baz2); 0 1 0; -sind(baz2) 0 cosd(baz2)];
tr_rotated = R * tr_rtv;
end


function [s_window] = get_s_window(tr_rotated, candidate_arr_windows, p_window_index)
n_candidate_arr_windows = size(candidate_arr_windows, 1);
candidate_s_window_energy = zeros(n_candidate_arr_windows - p_window_index,1);

tr_rotated_stack = sum(abs(tr_rotated(2:3,:)));

for i = (p_window_index + 1):n_candidate_arr_windows
    candidate_s_window_energy(i) = sum(tr_rotated_stack(candidate_arr_windows(i,1):candidate_arr_windows(i,2)).^2);
end

[~, s_window_index] = max(candidate_s_window_energy);
s_window = candidate_arr_windows(s_window_index,:);
end