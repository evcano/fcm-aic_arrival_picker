% Modified by Eduardo Valero Cano
% ---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------

function [t1] = compute_traveltime(z,v,sx,sz,rx,rz)
t1 = zeros(length(rz), 1);
t2 = t1;
for ii = 1:length(rz)
    % if source is deeper than receiver swap their positions
    if sz > rz(ii)
        sxtmp = rx;
        sztmp = rz(ii);
        rxtmp = sx;
        rztmp = sz;
    else
        sxtmp = sx;
        sztmp = sz;
        rxtmp = rx;
        rztmp = rz(ii);
    end
    [~,~,t1(ii)] = Current_1DRayTracing(z,v,sxtmp,sztmp,rxtmp,rztmp);
end
end