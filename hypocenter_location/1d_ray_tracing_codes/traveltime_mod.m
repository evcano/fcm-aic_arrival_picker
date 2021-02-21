% Modified by Eduardo Valero Cano
% ---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------

function [tp, ts] = traveltime_mod(vmod,sx,sz,rx,rz)
z = [0; vmod(:,1)];
vp = vmod(:,2);
vs = vmod(:,3);

tp = zeros(length(rz), 1);
ts = tp;

for r = 1:length(rz)
    % if source is deeper than receiver swap their positions
    if sz > rz(r)
        sxtmp = rx(r);
        sztmp = rz(r);
        rxtmp = sx;
        rztmp = sz;
    else
        sxtmp = sx;
        sztmp = sz;
        rxtmp = rx(r);
        rztmp = rz(r);
    end
    [~,~,tp(r)] = Current_1DRayTracing(z,vp,sxtmp,sztmp,rxtmp,rztmp);
    [~,~,ts(r)] = Current_1DRayTracing(z,vs,sxtmp,sztmp,rxtmp,rztmp);
end
end