function [x_min,f_min, iter] = BFGS(f, gf, x0, Stop_tol, c1, c2, varargin)
alpha_max = 100;
max_iter = 10000;
x_k = x0;
f_k = f(x0);
bracketing_eps = 1e-16;
% stop = false;
f_bar = 0;
delta_f = f_k - f_bar;

I = eye(length(x0));
c_k = I;
g_k = gf(x_k);
for iter = 1:max_iter
    p_k =  -c_k * g_k(:);
    % [alphaf, failed]  = LineSearch(phi_alpha, d_phi_alpha, alpha_max, c1, c2);
    
    Phi = @(alpha) f(x_k + alpha*p_k); 
    D_phi  = @(alpha) p_k'*gf(x_k + alpha*p_k);
    [alpha_k, phi_alpha_k, ~] = LineSearch(Phi, f_k , D_phi, p_k'*g_k, delta_f, c1, c2, alpha_max, bracketing_eps, varargin{1:end});
    x_k1 = x_k + alpha_k*p_k;
    
    g_k1 = gf(x_k1);
    f_k1 = phi_alpha_k;
    delta_f = f_k - f_k1;

    delta_k = x_k1 - x_k;
    gamma_k = g_k1 - g_k;
    g_k = g_k1;
    if iter == 1
        c_k = ((gamma_k'*delta_k)/(gamma_k'*gamma_k))*c_k;
    else
        dg_k = 1/(delta_k'*gamma_k);
        c_k = (I-dg_k*delta_k*gamma_k')*c_k*(I-dg_k*gamma_k*delta_k') + dg_k*(delta_k*delta_k');
    end
    
    x_k = x_k1;
    f_k = f_k1;
    
    n = norm(g_k1) ;
    if n <= Stop_tol
        break
    end
end
x_min = x_k;
f_min = f_k1;
end

