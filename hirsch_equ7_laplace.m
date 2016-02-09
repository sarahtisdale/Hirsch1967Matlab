function vprime = hirsch_equ7_laplace(R,v,params);
% Solve laplace (zero space charge) version of second order non-linear 
% differential equation 7 from Hirsch 1967 paper:
%   v''(R) + 2/R*v'(R) = 0
%
% Input:  
%   R: radius
%   v: potential (should be 'phi', but v is easier)
%
% Output:
%   2row 1col vector of v(r), v'(r)
% 
% Matlab ode45 only likes first derivatives, so let:  
%    v_1(R) = v(R)
%    v_2(R) = v'(R)
%
% We get this system of first order diff eqs:
%    v_1'(R) = v_2(R)
%    v_2'(R) = -2/R*v_2(R) 
%
% Get parameters K_plus and lambda_plus:
K_plus = params(1);
lambda_plus = params(2);

% build in pieces for easier debugging of problematic terms
vprime2 = -2/R*v(2);
vprime=[v(2); vprime2 ];
