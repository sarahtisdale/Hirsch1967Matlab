function vprime = hirsch_equ7(R,v,params);
% Solve second order non-linear differential equation 7 from Hirsch 1967 paper:
%   v''(R) + 2/R*v'(R) = K_plus/(R^2)*( v^(-1/2) - lamda_plus*(1-v)^(-1/2)
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
%    v_2'(R) = -2/R*v_2(R) + K_plus/(R^2)*( v_1(R)^(-1/2) - lamda_plus*(1-v_1(R))^(-1/2)
%
% Get parameters K_plus and lambda_plus:
K_plus = params(1);
lambda_plus = params(2);

% build in pieces for easier debugging of problematic terms
int1 = 1 - v(1);
int2 = int1^(-1/2);
int3 = v(1)^(-1/2);
int4 = K_plus/(R^2);
vprime2 = -2/R*v(2);
vprime2 = vprime2 + int4*( int3 - lambda_plus*(int2) );
vprime=[v(2); vprime2 ];
