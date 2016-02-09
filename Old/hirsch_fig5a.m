% Reproduce Figure 5 from 1967 Hirsch Inertial Electrostatic Confinement
% Fusion paper.
%
% NOTE: the Erratum makes important changes to this figure
%   - removes the dotted line series (lambda=0.3) plot
%   - clearly shows the plot is split in two pieces with different K+ &
%     lambda+ values
%

% clean slate
clc; clear all; close all;

% Convert R=r/r_c to R=r/r_a
%  - Figure 5 x-axis appears to be in terms of R=r/r_c
%  - But equ7 is written in terms of R=r/r_a
%  - Convert our x-axis values into to R=r/r_a before solving equ7
%  - Need to know the ratio of r_c/r_a.  I found them for K+=0.7 and
%    lambda+=0.454 in "hirsch_fig3_k07.m" 
% Actually, this R conversion seems to have no effect.  Plotting converted
% graph vs an unconverted one yields exactly the same graphs.  
%
% Now we use this solely as the location of the virtual anode to start/end 
% the different K/lambda portions of our graph.
Ra_over_Rc = 3.955;


% set up graph, axis labels and limits
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
legend_strs = cell(4,1);
pnum = 1;  % plot number 

%%%%%%%%%%%%
% outer K+=0.7 plot
%  - solved forward from virtual anode to real cathode
K_plus = 0.7;
lambda_plus = 0.454;
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

%%%%%%%%%%%%
% outer K+=0.7 plot
%  - solved backward from real cathode to virtual anode 
K_plus = 0.7;
lambda_plus = 0.454;
% radius varies from 1 down to virtual anode
R_c = 1:-0.005:(1/Ra_over_Rc);
% initial values: potential 1, slope zero
hinits=[0.999999999;0]; 
% solve
[R,v] = ode45(@hirsch_equ7,R_c,hinits,[],[K_plus lambda_plus]);
% plot
plot(R_c,v(:,1)); % only plot the first column (second column is v')
legend_strs{pnum} = sprintf('K_{+}=%04f, \\lambda_{+}=%0.3f, backward',...
                            K_plus, lambda_plus);
pnum = pnum+1;

%%%%%%%%%%%%
% inner K+=0.0859 plot
%  - solved forward from 0 to virtual anode
%  - ugly initial conditions
K_plus = 0.0859;
lambda_plus = 3.7;
% radius varies from virtual anode down to zero;
R_c = 0.01:0.005:1/Ra_over_Rc;
% initial values : messy!  potenial 1ish, slope negative infinity!?
hinits=[0.000001;99999]; 
% solve
[R,v] = ode45(@hirsch_equ7,R_c,hinits,[],[K_plus lambda_plus]);
% plot
plot(R_c,v(:,1)); % only plot the first column (second column is v')
legend_strs{pnum} = sprintf('K_{+}=%04f, \\lambda_{+}=%0.3f, forward',...
                             K_plus, lambda_plus);
pnum = pnum+1;


%%%%%%%%%%%%
% inner K+=0.0859 plot
%  - solved backward from virtual anode to zero
K_plus = 0.0859;
lambda_plus = 3.7;
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


% finish graph 
legend(legend_strs,'Location','southeast'); % add legend
hold off % no more plots for this figure


