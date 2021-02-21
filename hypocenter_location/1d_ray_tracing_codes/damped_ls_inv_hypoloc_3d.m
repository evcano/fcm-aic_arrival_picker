% ---------------------
% Author: Jubran Akram
% Supplementary material for our manuscript(A simple and robust
% implementation of two-point ray tracing for layered models with
% constant velocity) submitted for peer-review in Computers & Geosciences.
% February 7, 2016
%---------------------

function [x, res, flag,r_w, r] = damped_ls_inv_hypoloc_3d(vmod, r0, obT, p_nan, s_nan, W, x0, dx, dz, max_iter, stopping_criterion)
% This function computes the damped least square solution for hypocenter locations in three dimensional space
% dm = (G'G + A^2*I)*G'*r --- update (model, residual and G) -- stop when
% the stopping criterion is met or the maximum iterations are reached
%
%-------------- Input Parameters:
%
% 1. 1D velocity model = vmod = Isotropic case (num_layers x 3 [z, vp, vs])
% 2. Receiver positions = r0 = (num_observations x 3 [rec_x rec_y rec_z])
% 3. Observed arrival times = obT = P- and S-combined time vector
% 4. Missing arrival times = p_nan, s_nan = vectors
% 5. Weighting matrix = W = diagonal matrix
% 6. Initial guess = x0 = [x_0 y-0 z_0 t_0]'
% 7. Grid spacing for gradient computations = dx, dz = constants
% 8. Maximum number of iterations = max_iter = constant
% 9. Stopping criterion = stopping_criterion = constant
%
%-------------- Output Parameters:
%
% 1. Hypocenter location vector = x
% 2. Residual vector = res
% 3. Exit message = flag

% Parameter checks
if nargin < 7
    flag = strcat('You have entered ', num2str(nargin), ' parameters. Add the required 7 parameters');
    x = [];
    res = [];
    return
end

if nargin == 7
    dx = 1;
    dz = 1;
    max_iter = 25;
    stopping_criterion = 1e-9;
elseif nargin == 8
    dz = 1;
    max_iter = 25;
    stopping_criterion = 1e-9;
elseif nargin == 9
    max_iter = 25;
    stopping_criterion = 1e-9;
elseif nargin == 10
    stopping_criterion = 1e-9;
end

% main loop
ii = 1;
x = x0;

m = length(obT);
n = 4;

res = zeros(max_iter,1);
res(1) = 100000;
flag = "maximum iteration reached";

L = convmtx([-1 1],n-1);
D =   eye(n); %L'*L;
% num = 0;

p = length(nonzeros(p_nan));

M = ones(n,1);

res1 = 10000;
W1 = eye(m);
while ii <= max_iter
    
    res1_prev = res1;
    x_prev = x;
    
    r = compute_residuals(vmod, x, r0, obT, p_nan, s_nan);
    
%     if ii > 3
%     r1a =(abs(r(1:p)) - min(abs(r(1:p))))./(max(abs(r(1:p))) - min(abs(r(1:p))));
%     r1b =(abs(r(p+1:m-p)) - min(abs(r(p+1:m-p))))./(max(abs(r(p+1:m-p))) - min(abs(r(p+1:m-p))));
%     r1c =(abs(r(m-p+1:end)) - min(abs(r(m-p+1:end))))./(max(abs(r(m-p+1:end))) - min(abs(r(m-p+1:end))));
%     
%     r1 = [r1a; r1b; r1c] + eps;
%     W1 =  diag(r1.^2);
%     
%     end
    
    r_w = W1.*W*r;

    J = compute_jacobian(vmod,x,r0,dx, dz, p_nan, s_nan);
    
    J_w = W1.*W*J;

%     for kl = 1:n
%         M(kl) = norm(J_w(:,kl));
%         J_w(:,kl) = J_w(:,kl)./M(kl);
%         
%     end
    
   

    
    if ~isempty(nonzeros(isnan(J_w))); pn; end
    
    
    
    a = find_alpha(J_w,r_w,r,x,m,n,r0,p_nan,s_nan, obT,W1.*W,vmod,D);
    
    if isempty(a)
        
        
        flag = "residual cannot be decreased further";
        
        break;
        
    end
        
%     [a1 a]
    
    dx_sol = pinv(J_w'*J_w + a^2*D)*J_w'*r_w;
    
    x = x_prev + dx_sol;
    
    ii = ii + 1;
    
    res(ii) = sqrt(sum(r_w(1:m).^2)/(m-n));
    
    res1 = sqrt(sum(r(1:m).^2)/(m-n));    
    
    if abs(res1 - res1_prev) < 1e-7
        num = num +1;
        if num == 5
           flag = "residual couldn't be decreased any further"; 
           break;
        end
    else
        num = 0;
    end
    
    if res1 <= 0.001
       flag = "minimum residual criterion has met";
       break;
    end 
    
    
    
    fprintf("Iteration number %d  with a-value %f has the residual %f and time residual %f\n", ii-1, a, res(ii-1), res1);
    
end

end

function J = compute_jacobian(vmod,x,r0,dx, dz, p_nan, s_nan)
%% Computes Jacobian matrix for traveltimes in r-z plane
src_xp1 = x(1) + dx;
src_xm1 = x(1) - dx;
src_yp1 = x(2) + dx;
src_ym1 = x(2) - dx;
src_zp1 = x(3) + dz;
src_zm1 = x(3) - dz;

src_z = x(3);
rcv_x = 0.*r0(:,1);
rcv_z = r0(:,3);

azm_yp = atan2d((x(2)+dx-r0(1,2)),(x(1)-r0(1,1)));
azm_ym = atan2d((x(2)-dx-r0(1,2)),(x(1)-r0(1,1)));
azm_xp = atan2d((x(2)-r0(1,2)),(x(1)+dx-r0(1,1)));
azm_xm = atan2d((x(2)-r0(1,2)),(x(1)-dx-r0(1,1)));
azm_yp(azm_yp<0) = 360 + azm_yp(azm_yp<0);
azm_ym(azm_ym<0) = 360 + azm_ym(azm_ym<0);
azm_xp(azm_xp<0) = 360 + azm_xp(azm_xp<0);
azm_xm(azm_xm<0) = 360 + azm_xm(azm_xm<0);

% azm_yp = sind(azm_yp);
% azm_ym = sind(azm_ym);
% azm_xp = sind(azm_xp);
% azm_xm = sind(azm_xm);

ind_p = find(p_nan~=0);
ind_s = find(s_nan~=0);


p_length =length(ind_p);
s_length =length(ind_s);

num_recv = p_length + s_length;

J = ones(num_recv + p_length,4);

% P-S-rows
src_x = sqrt((src_xp1-r0(1,1)).^2 + (x(2)-r0(1,2)).^2);
[modP_xp1, modS_xp1] = traveltime_mod(vmod,src_x,src_z,rcv_x,rcv_z);
src_x = sqrt((src_xm1-r0(1,1)).^2 + (x(2)-r0(1,2)).^2);
[modP_xm1, modS_xm1] = traveltime_mod(vmod,src_x,src_z,rcv_x,rcv_z);

src_x = sqrt((x(1)-r0(1,1)).^2 + (src_yp1-r0(1,2)).^2);
[modP_yp1, modS_yp1] = traveltime_mod(vmod,src_x,src_z,rcv_x,rcv_z);
src_x = sqrt((x(1)-r0(1,1)).^2 + (src_ym1-r0(1,2)).^2);
[modP_ym1, modS_ym1] = traveltime_mod(vmod,src_x,src_z,rcv_x,rcv_z);

src_x = sqrt((x(1)-r0(1,1)).^2 + (x(2)-r0(1,2)).^2);
[modP_zp1, modS_zp1] = traveltime_mod(vmod,src_x,src_zp1,rcv_x,rcv_z);
[modP_zm1, modS_zm1] = traveltime_mod(vmod,src_x,src_zm1,rcv_x,rcv_z);

J(1:p_length,1) = (modP_xp1(ind_p) - modP_xm1(ind_p))/(2*dx);
J(1:p_length,2) = (modP_yp1(ind_p) - modP_ym1(ind_p))/(2*dx);
J(1:p_length,3) = (modP_zp1(ind_p) - modP_zm1(ind_p))/(2*dz);

J(p_length+1:num_recv,1) = (modS_xp1(ind_s) - modS_xm1(ind_s))/(2*dx);
J(p_length+1:num_recv,2) = (modS_yp1(ind_s) - modS_ym1(ind_s))/(2*dx);
J(p_length+1:num_recv,3) = (modS_zp1(ind_s) - modS_zm1(ind_s))/(2*dz);

J(num_recv+1:num_recv+p_length,1) = (azm_xp - azm_xm)/(2*dx);
J(num_recv+1:num_recv+p_length,2) = (azm_yp - azm_ym)/(2*dx);
J(num_recv+1:num_recv+p_length,[3 4]) = 0;

end

function r = compute_residuals(vmod, x, r0, obT, p_nan, s_nan)
%% This function computes traveltime residuals (r)
src_x = sqrt((x(1)-r0(1,1)).^2 + (x(2)-r0(1,2)).^2);
src_z = x(3);
rcv_x = 0.*r0(:,1);
rcv_z = r0(:,3);
azm = atan2d((x(2)-r0(1,2)),(x(1)-r0(1,1)));
azm(azm<0) = 360 + azm(azm<0);
% azm = sind(azm);
[modP, modS] = traveltime_mod(vmod,src_x,src_z,rcv_x,rcv_z);
modT = [modP(p_nan ~=0)+x(4); modS(s_nan ~=0)+x(4); azm.*ones(length(nonzeros(p_nan)),1)];

r = obT-modT;
end

function a = find_alpha(J,r,r1, x, m, n, r0,p_nan, s_nan, obT,W,vmod,D)
%% This function finds the regularization parameter
[U,S,V] = svd(J);
s_i = diag(S);
a_length = 100;

a_range = sort(log10(0.001) + (log10(5000)-log10(0.001))*rand(a_length,1));
% a_range = linspace(log10(min(s_i)), log10(max(s_i)), a_length);

a_range = 10.^(a_range);

% L = zeros(a_length,1);
% m_new = L;
% P = L;
p = length(nonzeros(p_nan));
res = sqrt(sum(r(1:m).^2)/(m-n));
q = zeros(a_length,1);

parfor ii = 1:a_length
    
    J_plus = pinv(J'*J + a_range(ii)^2*D)*J';
    m_plus = J_plus*r;
    
    r_new = compute_residuals(vmod, x + m_plus, r0, obT, p_nan, s_nan);
        
    r_w = W*r_new;
    
    q(ii) = sqrt(sum(r_w(1:m).^2)/(m-n));
    
%     if sqrt(sum(r_w.^2)/(m-n)) < res
%         a = a_range(ii);
%         break;
%         
end
ind = find(res - q > 0, 1, 'first');

a = a_range(ind);





end


