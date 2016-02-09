% Charge density for K_plus = 1.0 obtained by first solving 
% equation 7 for v(R).  SI Units were used by incorporating 
% the elctric constant K_E = 1 / (4*pi*eps_0) into the definition
% of K_plus, to give K_plus = K_E*4*pi*r^2*rho*v^1/2 / |V_o|
% solving for the charge density gives:
% rho = eps_0*|V_o|/(r^2*v^1/2)

% choose parameters
K_plus = 1.0;
lambda_plus = 0;

% radius from 1 to 4
R = 1:0.01:4;

% intial values at radius 1 (from the text):
%   v  = 0.0000001  - nudge slightly above zero to avoid infinity NaNs
%   v' = 0
hinits=[0.0000001,0];

% solve for v
[R,v] = ode45(@hirsch_equ7,R,hinits,[],[K_plus lambda_plus]);

% let r = 0 correspond to the location of the anode
% and r = 1 correspond to the location of the cathode
r = linspace(0,1,length(R));  % distance from the virtual anode (m)
r = r'; % transpose r into a column vector
eps_0 = 8.85*10^-12;  % permittivity of free space (F/m)
V_o = -1000;  % applied voltage (V)
rho = eps_0*abs(V_o)/(r.^2 .* v(:,1).^(1/2));  % volume charge density (C/m^3)

% axis labels and limits
figure(2)
plot(r,rho)
title('Charge Density between Virtual Anode and Cathode');
ylabel('\rho(r)');
xlabel('r');
ax = gca;
ax.XTick = 0:.1:1;
ax.XLim = [0 1];
ax.YTick = 0:0.2E-9:2E-9;
ax.YLim = [0 2E-9];
ax.PlotBoxAspectRatio = [1,0.635,1];







