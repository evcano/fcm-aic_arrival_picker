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

function event_results = verify_phase_labeling(event_results, par)
n_receivers = length(event_results.p_pick);
rno(:,1) = 1:n_receivers;

% RANSAC parameters
sampleSize = par.ransac_samplesize;
maxDistance = par.ransac_maxdistance;
fitFunction = @(data) polyfit(data(:,1),data(:,2),2);
misfitFunction = @(model, data) sum((data(:,2) - polyval(model,data(:,1))).^2,2);

% unknown-arrival event
if sum(event_results.pick_flag == 0) == n_receivers
    return;
end

% fit s-wave moveout
u_rec = find(event_results.pick_flag == 0);
x(:,1) = event_results.s_pick;
x(u_rec) = event_results.u_pick(u_rec);
tmp = isnan(x);

validateFunction = @(model) verify_s_moveout(model,rno(~tmp),event_results.p_pick(~tmp));
[modelRANSAC, inliers] = ransac([rno(~tmp),x(~tmp)], fitFunction, misfitFunction,...
    sampleSize, maxDistance, 'ValidateModelFcn',validateFunction);
event_results.s_moveout = round(polyval(modelRANSAC,rno));
event_results.s_moveout_ni = sum(inliers);

% label unknown arrivals
for i = 1:length(u_rec)
    r = u_rec(i);
    dif = event_results.s_moveout(r) - event_results.u_pick(r);
    if abs(dif) <= maxDistance || dif < 0
        event_results.s_pick(r) = event_results.u_pick(r);
        event_results.s_window(r,:) = event_results.u_window(r,:);
        event_results.u_pick(r) = NaN;
        event_results.u_window(r,:) = NaN(length(r),2);
        event_results.pick_flag(r) = 2;
    else
        event_results.p_pick(r) = event_results.u_pick(r);
        event_results.p_window(r,:) = event_results.u_window(r,:);
        event_results.u_pick(r) = NaN;
        event_results.u_window(r,:) = NaN(length(r),2);
        event_results.pick_flag(r) = 1;
    end
end

% correct s picks identified as p picks
dif = abs(event_results.s_moveout - event_results.p_pick);
tmp = find(dif < maxDistance);
event_results.s_pick(tmp) = event_results.p_pick(tmp);
event_results.s_window(tmp,:) = event_results.p_window(tmp,:);
event_results.p_pick(tmp) = NaN;
event_results.p_window(tmp,:) = NaN(length(tmp),2);
event_results.pick_flag(tmp) = 3;

% fit p-wave moveout
tmp = isnan(event_results.p_pick);
if sum(tmp) < 8
    x(:,1) = event_results.p_pick;
    
    validateFunction = @(model) verify_p_moveout(model,rno(~tmp),event_results.s_pick(~tmp));
    [modelRANSAC,inliers] = ransac([rno(~tmp),x(~tmp)], fitFunction, misfitFunction,...
        sampleSize, maxDistance, 'ValidateModelFcn',validateFunction);
    
    event_results.p_moveout = round(polyval(modelRANSAC,rno));
    event_results.p_moveout_ni = sum(inliers);
end
end

function isValid = verify_p_moveout(model,varargin)
y = varargin{1};
s_pick = varargin{2};
p_moveout = round(polyval(model,y));
dif = s_pick - p_moveout;
if all(dif(~isnan(dif)) >= 0)
    isValid = true;
else
    isValid = false;
end
end

function isValid = verify_s_moveout(model,varargin)
y = varargin{1};
p_pick = varargin{2};
s_moveout = round(polyval(model,y));
dif = s_moveout - p_pick;
if all(dif(~isnan(dif)) >= 0)
    isValid = true;
else
    isValid = false;
end
end