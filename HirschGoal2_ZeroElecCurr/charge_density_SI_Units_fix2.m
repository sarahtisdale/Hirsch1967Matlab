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
%R = 1:0.01:4;  %sst: wrong radius range.  should be 0 to 1



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% new figure for fig 5

% Location of the virtual anode to start/end the different K/lambda 
% portions of our graph. (see hirsch_fig3_k10.m for how this is calculated)
Ra_over_Rc = 2.315;

% set up graph, axis labels and limits
figure('NumberTitle','off','Name','Figure 5');
hold on % several graphs on same figure
title('Fig 5. Normalized potential distribution for K_{+}=0.7 and \lambda_{+}=lambda_{+max}=0.454 at the real cathode.');
ylabel('\phi(R)');
xlabel('R');
ax = gca;
ax.XTick = 0:0.1:1;
ax.XLim = [0 1];
ax.YTick = 0:0.2:1;
ax.YLim = [-.05 1];
ax.PlotBoxAspectRatio = [1,0.51,1];  % copy aspect ratio of the paper's graph
grid on
legend_strs = cell(2,1);
pnum = 1;  % plot number 

%%%%%%%%%%%%
% outer K+=1.0 plot
%  - solved forward from virtual anode to real cathode
K_plus = 1.0;
lambda_plus = 0;
% radius varies from virtual anode to 1
R_c = 1/Ra_over_Rc:0.005:1;
% initial values: potential and slope are zero
hinits=[0.00000000001;0]; 
% solve
[R,v] = ode45(@hirsch_equ7,R_c,hinits,[],[K_plus lambda_plus]);
% plot
plot(R_c,v(:,1)); % only plot the first column (second column is v')
legend_strs{pnum} = sprintf('K_{+}=%04f, \\lambda_{+}=%0.3f, forward', ...
                             K_plus, lambda_plus);
pnum = pnum+1;
Rsave1 = R;
Vsave1 = v;

%%%%%%%%%%%%
% inner K+=0.0859 plot
%  - PROBLEM!!!!!  
%      Don't know what "inner K+/lambda+" to use when outer K+=1.0 
%      and outer lambda+=0.  My guesses:  K+=0.1, lambda+=0
%            
%  - solved backward from virtual anode to zero
K_plus = 0.1;
lambda_plus = 0;
% radius varies from virtual anode down to zero;
R_c = 1/Ra_over_Rc:-0.005:0;
% initial values : potenial and slope zero
hinits=[0.0000001;0]; 
% solve
[R,v] = ode45(@hirsch_equ7,R_c,hinits,[],[K_plus lambda_plus]);
% plot
plot(R_c,v(:,1)); % only plot the first column (second column is v')
legend_strs{pnum} = sprintf('K_{+}=%04f, \\lambda_{+}=%0.3f, backward',...
                             K_plus, lambda_plus);
pnum = pnum+1;
Rsave2=R;
Vsave2=v;

% finish graph 
legend(legend_strs,'Location','southeast'); % add legend
hold off % no more plots for this figure

% set up full R&V from two previous plots
r = [flipud(Rsave2); Rsave1];
v = [flipud(Vsave2); Vsave1];


%%% josh's stuff below here vvvvv (edited lines noted)

% let r = 0 correspond to the location of the anode
% and r = 1 correspond to the location of the cathode

%sst: use r and v from above fig 5 data)
%r = linspace(0,1,length(R));  % distance from the virtual anode (m)
%r = r'; % transpose r into a column vector

eps_0 = 8.85*10^-12;  % permittivity of free space (F/m)
V_o = -1000;  % applied voltage (V)

%sst: add .* and ./ to fix rho
rho = eps_0.*abs(V_o)./(r.^2 .* v(:,1).^(1/2));  % volume charge density (C/m^3)

% axis labels and limits
figure(2)
plot(r,rho)
title('Charge Density between Virtual Anode and Cathode');
ylabel('\rho(r)');
xlabel('r');
ax = gca;
ax.XTick = 0:.1:1;
ax.XLim = [0 1];

%sst: fix Y limits, Y scale, and aspect ratio
ax.YLim = [10^-9 10^-3];
ax.YScale = 'log';
ax.PlotBoxAspectRatio = [1,1.2,1];
grid on;







