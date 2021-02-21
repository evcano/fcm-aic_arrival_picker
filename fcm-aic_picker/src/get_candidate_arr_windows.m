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

function [candidate_arr_windows] = get_candidate_arr_windows(signal_membership, signal_threshold, tdom)
signal_membership(signal_membership < signal_threshold) = 0;
signal_membership(signal_membership >= signal_threshold) = 1;

limits = diff(signal_membership);

if signal_membership(1) ~= 0
    limits = [1 limits];
else
    limits = [0 limits];
end

if signal_membership(end) ~= 0
    if signal_membership(end -1) == 1
        limits(end) = -1;
    else
        limits(end) = 0;
    end
end

candidate_arr_windows(:,1) = find(limits == 1);
candidate_arr_windows(2:end,1) = candidate_arr_windows(2:end,1) -1;
candidate_arr_windows(:,2) = find(limits == -1);

candidate_arr_windows((candidate_arr_windows(:,2)-candidate_arr_windows(:,1)) < tdom,:) = [];
end