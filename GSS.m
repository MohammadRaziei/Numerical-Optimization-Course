function [x_min, N] = GSS(f, a, b, epsilon, varargin)
% Golden Selection Search 
% created by @MohammadRaziei
% 
% f    : target function
% a, b : intervals
% b    : tolerance
% N    : number of iteration
% % % % % % % % % % % % % % % % % % % % % % % %
if(a > b), c = a; a = b; b = c; clear c; end
if(nargin < 4)
    epsilon = 1e-8;
end

rho = (3 - sqrt(5))/2;

x1 = a + rho * (b-a);
x2 = a + (1-rho) * (b-a);

N = 0;
while true
    N = N + 1;
    % calc f(a), f(b)
    f1 = f(x1, varargin{1:end});
    f2 = f(x2, varargin{1:end});
    % define intervals
    if(f1 < f2)
        b = x2;
        x2 = x1;
        x1 = a + rho * (b-a);
    else
        a = x1;
        x1 = x2;
        x2 = a + (1-rho) * (b-a);
    end
    
    if(norm(x2 - x1) < epsilon)
        break
    end
    
end

x_min = mean([x1 x2]);

% Mohammad Raziei
% 98206223
% Thursday, October 1, 2020

end