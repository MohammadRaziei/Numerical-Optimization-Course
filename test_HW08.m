clc; clear all; close all;
% By @MohammadRaziei
%% Rosenbrock 
clc; f_rosenbrock; gf_rosenbrock; % reset counters
stop_tol = 1e-5; c1 = 1e-4; c2 = 0.9; 
f = @(x) f_rosenbrock(x(1), x(2));
gf = @(x) gf_rosenbrock(x(1), x(2));
x0 = [1;2];
[rosenbrock_min_x, rosenbrock_min_f, rosenbrock_min_iter] = BFGS(f, gf, x0, stop_tol,c1,c2);
rosenbrock_counter_f = f_rosenbrock();
rosenbrock_counter_gf = gf_rosenbrock();

%% Powel
clc; f_powel; gf_powel; % reset counters
stop_tol = 1e-5; c1 = 1e-4; c2 = 0.9; 
f = @(x) f_powel(x(1), x(2), x(3), x(4));
gf = @(x) gf_powel(x(1), x(2), x(3), x(4));
x0 = [1;2;2;2];
[powel_min_x, powel_min_f, powel_min_iter] = BFGS(f, gf, x0, stop_tol,c1,c2);
powel_counter_f = f_powel();
powel_counter_gf = gf_powel();

%% Create report
Title = "HW08 (BFGS Algorithm)"; filename = 'README_HW08.md';
s = @(a)"["+strjoin(string(a(:)),', ')+"]";
readme = [
"# "+Title+"  ";
"  ";
"  ";
"## BFGS  ";
"  ";
"|   | final x  | final f | # func iter | # func eval | # grad eval | # Hess eval  |";
"|:---:|:---:|:---:|:---:|:---:|:---:|:---:|";
"| Powel | "+s(powel_min_x)+" | "+powel_min_f+" | "+powel_min_iter+" | "+powel_counter_f+" | "+powel_counter_gf+" | - |";
"| Rosenbrock | "+s(rosenbrock_min_x)+" | "+rosenbrock_min_f+" | "+rosenbrock_min_iter+" | "+rosenbrock_counter_f+" | "+rosenbrock_counter_gf+" | - |";
"  ";];
fileID = fopen(filename,'w'); [~]= fprintf(fileID,'%s\n',readme); [~]= fclose(fileID);
clear fileID s filename Title;
winopen('README_HW08.md');
%% FUNCTIONS
% %% CALC Rosenbrock
% clc; syms x1 x2
% f = 100*(x2-x1^2)^2 + (1-x1)^2
% gf = [diff(f, x1);diff(f, x2)]
% Hf = [diff(gf,x1), diff(gf,x2)]
% %% CALC Powel
% clc; syms x1 x2 x3 x4
% f = (x1 + 10*x2)^2 + 5*(x3 - x4)^2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4;
% gf = [diff(f, x1); diff(f, x2); diff(f, x3); diff(f, x4)]
% Hf = [diff(gf,x1), diff(gf,x2), diff(gf,x3), diff(gf,x4)]
clear f gf
%% Rosenbrock function

function y = f_rosenbrock(x1, x2)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = 100*(x2-x1^2)^2 + (1-x1)^2;
end
end

function y = gf_rosenbrock(x1, x2)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = [2*x1 - 400*x1*(- x1^2 + x2) - 2;   - 200*x1^2 + 200*x2];
%     x = [x1 x2];
%     y = [-400*(x(2) - x(1)^2)*x(1) - 2*(1 - x(1)) ; 200*(x(2) - x(1)^2)];
end
end

function y = Hf_rosenbrock(x1, x2)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = [ 1200*x1^2 - 400*x2 + 2, -400*x1;  -400*x1,     200];
end
end



%% Powel function

function y = f_powel(x1, x2, x3, x4)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = (x1 + 10*x2)^2 + 5*(x3 - x4)^2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4;
end
end

function y = gf_powel(x1, x2, x3, x4)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = [2*x1 + 20*x2 + 40*(x1 - x4)^3;
        20*x1 + 200*x2 + 4*(x2 - 2*x3)^3;
        10*x3 - 10*x4 - 8*(x2 - 2*x3)^3;
        10*x4 - 10*x3 - 40*(x1 - x4)^3];
end
end

function y = Hf_powel(x1, x2, x3, x4)
persistent counter; if isempty(counter), counter = 0; end
if and(nargout == 0, nargin == 0), counter = 0; return ;end %% reset counter
if(nargin == 0), y = counter; else, counter = counter + 1;
    y = [ 120*(x1 - x4)^2 + 2,                     20,                     0,     -120*(x1 - x4)^2;
                           20, 12*(x2 - 2*x3)^2 + 200,     -24*(x2 - 2*x3)^2,                    0;
                            0,      -24*(x2 - 2*x3)^2, 48*(x2 - 2*x3)^2 + 10,                  -10;
             -120*(x1 - x4)^2,                      0,                   -10, 120*(x1 - x4)^2 + 10];
end
end