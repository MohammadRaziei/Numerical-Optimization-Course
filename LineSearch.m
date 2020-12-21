function [alpha, phi_alpha, stop] = LineSearch(Phi, phi_0, D_phi, d_phi_0, delta_f, c1, c2, alpha_max, bracketing_eps, varargin)
% function [alpha] = LS_SW(Phi , D_phi , p , phi_0 , alpha , alpha_1)
%     if nargin < 8, zoom_eps = 1e-16; end
%     if nargin < 9, gss_tol = 1e-5; end
   
%     c1 = 10^-3;
%     c2 = 0.5;1
    stop = false;
    alpha = 0;
    alpha_0 = 0;
%     bracketing_eps = 1e-16;
    i = 1;
    alpha_i_1 = alpha_0;
    alpha_i = ChooseAlpha1(delta_f, d_phi_0);
    d_phi_alpha_i_1 = d_phi_0;
    phi_alpha_i_1 = phi_0;

    while true
%         alpha_i
        phi_alpha_i = Phi(alpha_i);
        if (phi_alpha_i > phi_0 + c1*alpha_i* d_phi_0) || ((i>1) && phi_alpha_i >= phi_alpha_i_1)
            [alpha, phi_alpha] = ZoomSelection(alpha_i_1 , alpha_i , Phi , D_phi, c1, c2, d_phi_alpha_i_1 , phi_0 , d_phi_0 , phi_alpha_i_1 , phi_alpha_i);
            break;
        end
        d_phi_alpha_i = D_phi(alpha_i);
        if abs(d_phi_alpha_i) < bracketing_eps
            stop = true;
            alpha = alpha_i;
            phi_alpha = phi_alpha_i;
%             d_phi_alpha = d_phi_alpha_i;
            break;
        end
        if (abs(d_phi_alpha_i) <= -c2 * d_phi_0)
            alpha = alpha_i;
            phi_alpha = phi_alpha_i; 
%             d_phi_alpha = d_phi_alpha_i;
            break;
        end
        if (d_phi_alpha_i >= 0)
%             [alpha, phi_alpha, d_phi_alpha] = ZoomSelection(alpha_i_1 , alpha_i , Phi , D_phi, c1, c2, d_phi_alpha_i_1 , phi_0 , d_phi_0 , phi_alpha_i_1 , phi_alpha_i);
            [alpha, phi_alpha] = ZoomSelection(alpha_i , alpha_i_1 , Phi , D_phi, c1, c2, d_phi_alpha_i_1 , phi_0 , d_phi_0 , phi_alpha_i_1 , phi_alpha_i);
            break;
        end
        i = i + 1;
        alpha_i_new = ChooseAlphaNew(alpha_i, alpha_i_1, alpha_max);
        
        alpha_i_1 = alpha_i;
        alpha_i = alpha_i_new;
        d_phi_alpha_i_1 = d_phi_alpha_i;
        phi_alpha_i_1 = phi_alpha_i;
    end

end
    
function alpha_1 = ChooseAlpha1(delta_f, d_phi_0)
    alpha_1 =  min(1,abs(-2.02*delta_f/d_phi_0));
end

function alpha = ChooseAlphaNew(alpha_i, alpha_i_1, alpha_max)
    T1 = 10;
    if alpha_max <= (2*alpha_i - alpha_i_1)
        alpha = alpha_max;
    else
        alpha = mean([(2*alpha_i - alpha_i_1) ,  min(alpha_max, alpha_i+T1*(alpha_i - alpha_i_1))]);
    end
end

function [alpha, phi_alpha] = ZoomSelection(alpha_i_1 , alpha_i , Phi , D_Phi, c1, c2, d_phi_alpha_i_1 , phi_0 , d_phi_0 , phi_alpha_i_1 , phi_alpha_i)

    zoom_eps = 1e-16;
    while true
        candid_a = alpha_i_1;   phi_condid_a = phi_alpha_i_1;
        candid_b = alpha_i;     phi_condid_b = phi_alpha_i;
        d_phi_candid_a = d_phi_alpha_i_1;
        condid_alpha = ChooseAlphaInZoom(candid_a, candid_b, phi_condid_a , phi_condid_b , d_phi_candid_a);
        phi_condid_alpha = Phi(condid_alpha); 

%         d_phi_condid_alpha = D_Phi(condid_alpha);
        if (candid_a - candid_b) * d_phi_candid_a < zoom_eps
            alpha = condid_alpha;
            phi_alpha = phi_condid_alpha; 
%             d_phi_alpha = d_phi_condid_alpha;
            break;
        end
        if (phi_condid_alpha > phi_0 + c1*condid_alpha*d_phi_0) || (phi_condid_alpha >= phi_alpha_i_1)
            candid_b = condid_alpha;
        else
            d_phi_condid_alpha = D_Phi(condid_alpha); 
            if (abs(d_phi_condid_alpha) <= -c2*d_phi_0)
                alpha = condid_alpha;
                phi_alpha = phi_condid_alpha;  
%                 d_phi_alpha = d_phi_condid_alpha;
                break;
            end
            candid_a_new = candid_a;
            candid_a = condid_alpha;
            d_phi_candid_a = d_phi_condid_alpha;
            if ((candid_b-candid_a)*d_phi_condid_alpha >= 0)
                candid_b = candid_a_new;
%             else
                %candid_b = candid_b;
            end
        end
        
        
        
        if abs(candid_a-candid_b)/abs(alpha_i_1-alpha_i) < 0.5
            alpha = condid_alpha;
            phi_alpha = phi_condid_alpha;  
%             d_phi_alpha = d_phi_condid_alpha;
            alpha_i_1 = candid_a;
            d_phi_alpha_i_1 = d_phi_candid_a;
            alpha_i = candid_b;
        else  
            alpha_m = (alpha_i_1+alpha_i)/2;
            phi_alpha_m = Phi(alpha_m);
%             d_phi_alpha_m = D_Phi(alpha_m);
            if ((alpha_i_1 - alpha_i) * d_phi_alpha_i_1 < zoom_eps)
                alpha = alpha_m;
                phi_alpha = phi_alpha_m; 
%                 d_phi_alpha = d_phi_alpha_m;
                break;
            end
            if (phi_alpha_m > phi_0 + c1*alpha_m*d_phi_0) || (phi_alpha_m >= phi_alpha_i_1)
                alpha_i = alpha_m;
            else
                d_phi_alpha_m = D_Phi(alpha_m);
                if (abs(d_phi_alpha_m) <= -c2*d_phi_0)
                    alpha = alpha_m;
                    phi_alpha = phi_alpha_m; 
%                     d_phi_alpha = d_phi_alpha_m;
                    break;
                end
                a_new = alpha_i_1;
                alpha_i_1 = alpha_m;
                d_phi_alpha_i_1 = d_phi_alpha_m;
                if (alpha_i - alpha_i_1) * d_phi_alpha_m >= 0
                    alpha_i = a_new;
%                 else
                    %b = b;
                end
            end
        end
    end

end




function  alpha_j = ChooseAlphaInZoom(a,b, phi_a, phi_b, d_phi_a)
%     gss_tol = 1e-5;
    phi0 = phi_a;
    phi1 = phi_b;
    d_phi0 = (b-a)*d_phi_a;

    z_min = -d_phi0 / (2*(phi1-phi0-d_phi0));
    alpha_j = a + z_min*(b-a);
    if (alpha_j > a && alpha_j > b) || (alpha_j < a && alpha_j < b)
        alpha_j = (a+b)/2;
    end
end
