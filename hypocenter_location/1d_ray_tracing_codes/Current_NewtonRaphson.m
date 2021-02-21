% ---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------

function q = Current_NewtonRaphson(h,hm,ek,delX,q0, niter)
%---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------
np = ((ek).^2.*(1+(hm/q0)^2) - 1).^(1/2);
f = (h'*(1./np) - delX);
f1 = h'*((ek).^2.*(hm^2/q0^3)./np.^3);
delQ = -f./f1;

i = 0;

while abs(f) > 1E-2
    
    q = q0 + delQ;
    
    np = ((ek).^2.*(1+(hm/q)^2) - 1).^(1/2);
    f = (h'*(1./np) - delX);
    f1 = h'*((ek).^2.*(hm^2/q^3)./np.^3);
    
    delQ = -f./f1;
    q0 = q;
    
    i = i + 1;
    
    if i== niter
        
        break;
        
    end
    
end

q = h.*(1./np);
