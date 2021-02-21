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

function [U_new, V, it] = fuzzy_c_means(X, n_clusters, n_iterations, fuzzifier, stop_criteria)
e = -2 / (fuzzifier - 1); % exponent
[n_data, n_features] = size(X);

V_ini = rand(n_clusters,n_features); % initialize centroids randomly
V = zeros(n_clusters,n_features);
U = zeros(n_data,n_clusters);
U_new = zeros(n_data,n_clusters);

% initialize membership matrix
b = 0;
for j = 1:n_clusters
    b = b + sqrt(sum((X - repmat(V_ini(j,:),n_data,1)).^ 2, 2)).^e;
end
for i = 1:n_clusters
    a = sqrt(sum((X - V_ini(i,:)).^2, 2)).^e;
    U(:,i) = a ./ b;
end

for it = 1:n_iterations
    % update centroids
    for i = 1:n_clusters
        u = U(:,i).^fuzzifier;
        V(i,:) = sum(u .* X) / sum(u);
    end
    
    % update membership matrix
    b = 0;
    for j = 1:n_clusters
        b = b + sqrt(sum((X - V(j,:)).^2, 2)).^e;
    end
    for i = 1:n_clusters
        a = sqrt(sum((X - V(i,:)).^2, 2)).^e;
        U_new(:,i) = a ./ b;
    end
    
    % evaluate misfit function
    misfit = 0;
    for i = 1:n_clusters
        misfit = misfit + (U_new(:,i).^fuzzifier) .* (sum((X - V(i,:)).^2, 2));
    end
    misfit = sum(misfit);
    
    % evaluate stop criterion
    dif = norm(U_new - U);
    if dif < stop_criteria
        break
    else
        U = U_new;
    end
end
end