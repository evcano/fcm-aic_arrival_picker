%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.


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