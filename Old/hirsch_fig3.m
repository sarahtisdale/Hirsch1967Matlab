% Reproduce Figure 3 from 1967 Hirsch Inertial Electrostatic Confinement
% Fusion paper.

% clean slate
clc; clear all; close all;

% radius from 1 to 4
R = 1:0.01:4;

% intial values at radius 1 (from the text):
%   v  = 0.0000001  - nudge slightly above zero to avoid infinity NaNs
%   v' = 0
hinits=[0.0000001,0];

% K_plus constant 1.0 for figure 3
K_plus = 1.0;
% lambda_plus values
lambda_plusses = [0,0.1,0.3,0.5,0.55,0.6,0.8];
legend_strs = cell(numel(lambda_plusses),1);

% several lambda_plus graphs on same graph
hold on

% axis labels and limits
title('Fig 3. Potential distributions at K_{+}=1.0.');
ylabel('\phi(R)');
xlabel('R');
ax = gca;
ax.XTick = 0:1:4;
ax.XLim = [0 4];
ax.YTick = 0:0.2:1;
ax.YLim = [0 1];
ax.PlotBoxAspectRatio = [1,0.635,1];

grid on

% solving & plotting
label_idx = round((numel(R)-1)/3)+1;
for i=1:numel(lambda_plusses),
    lambda_plus = lambda_plusses(i);
    legend_strs{i} = sprintf('\\lambda_{+}=%0.2f',lambda_plus);
    [R,v] = ode45(@hirsch_equ7,R,hinits,[],[K_plus lambda_plus]);
    plot(R,v(:,1)); % only plot the first column (second column is v')

    % create a silly text label
    alignment='left';
    label_str = sprintf('\\leftarrow %s',legend_strs{i});
    if logical(mod(i,2))
        alignment='right';
        label_str = sprintf('%s \\rightarrow',legend_strs{i});
    end
    text(R(label_idx),v(label_idx),label_str,'HorizontalAlignment',alignment);
end

% add a legend
legend(legend_strs,'Location','southeast');

% no more plots for this figure
hold off

