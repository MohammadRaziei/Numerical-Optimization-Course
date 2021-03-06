function [x_min, f_min, iter] = Newton_GSS(f, gf, Hf, x0, stop_tol, gss_tol, varargin)
a = 0; b = 10;
max_iter = 10000;
x_k = x0;
for iter = 1 : max_iter
    H_k = Hf(x_k);
    g_k = gf(x_k);
    p_k = H_k\(-g_k);
    
    alpha_k = GSS(@(alpha) f(x_k + alpha*p_k), a, b, gss_tol, varargin{1:end});
    x_k1 = x_k + alpha_k * p_k;
    if norm(x_k - x_k1) < stop_tol
        x_k = x_k1;
        break
    end
    x_k = x_k1;
end

x_min = x_k;
f_min = f(x_min);
end