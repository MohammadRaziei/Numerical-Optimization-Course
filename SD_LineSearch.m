function [x_min, f_min, iter] = SD_LineSearch(f, gf, x0, stop_tol, c1, c2, varargin)
    if nargin <  4, stop_tol = 1e-3;  end
    if nargin <  5, c1 = 1e-5;        end
    if nargin <  6, c2 = 1e-5;        end
    alpha_max = 100;

    max_iter = 10000;
    x_k = x0;
    f_k = f(x0);
    bracketing_eps = 1e-16;
    stop = false;
    f_bar = 0;
    delta_f = f_k - f_bar;
    for iter = 1 : max_iter
        g_k = gf(x_k);
%         if stop && all(abs(g_k) < bracketing_eps), break; end
        p_k = -g_k(:);

        

        Phi = @(alpha) f(x_k + alpha*p_k);
        D_phi  = @(alpha) p_k'*gf(x_k + alpha*p_k);
        [alpha_k, phi_alpha_k, stop] = LineSearch(Phi, f_k , D_phi, p_k'*g_k, delta_f, c1, c2, alpha_max, bracketing_eps, varargin{1:end});

        %         if alpha_k == 0
        x_k1 = x_k + alpha_k * p_k;
        f_k1 = phi_alpha_k;
%         g_k = d_phi_alpha_k;
        
        delta_f = f_k - f_k1;

        if (norm(x_k - x_k1) < stop_tol ) || stop
            x_k = x_k1;
            f_k = f_k1;
            break
        end
        x_k = x_k1;
        f_k = f_k1;
    end

    x_min = x_k;
    f_min = f_k;
end