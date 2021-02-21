%---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------

function [x,z,t] = Current_1DRayTracing(zmod,vmod,sx,sz,rx,rz)
n = length(zmod);

delX = rx - sx;

%--------Layering index ----------
indz = 1:n;

ind1 = interp1(zmod,indz,sz);
ind2 = interp1(zmod,indz,rz);

sind = floor(ind1);
rind = floor(ind2);

%--------Computing raypaths and time ---------

if sind == rind
    
    v = vmod(sind);
    
    z = [sz rz];
    x = [sx rx];
    t = sqrt((sz-rz)^2+(sx-rx)^2)/v;
    
else
    
    
    if zmod(rind) == rz
        
        rind = max(1,rind-1);
        
    end
    
    v = vmod(sind:rind);
    
    z = vertcat(sz,zmod(sind+1:rind),rz);
    
    
    
    n = length(z);
    
    [vm,M] = max(v);
    
    ek = vm./v;
    
    h = diff(z);
    
    a = h'*(1./(h(M)*(ek)));
    b = h(1:M-1)'*(1./sqrt((ek(1:M-1).^2)-1))+...
        h(M+1:end)'*(1./sqrt((ek(M+1:end).^2)-1));
    
    if isempty(b)
        b = a;
    end
    
    delC = a*b/(a-1);
    
    nx = abs(delX);
    
    if nx <= delC
        
        q0 = nx/a;
        
    else
        
        q0 = nx - b;
        
    end
    
    
    
    q = Current_NewtonRaphson(h,h(M),ek,nx,q0, 50) ;
    
    t = sqrt(q.^2+h.^2)'*(1./v);
    
    q = [0;q];
    
    x = sx + sign(delX)*tril(ones(n,n))*q;
    
end
