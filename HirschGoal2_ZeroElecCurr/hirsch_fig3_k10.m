% Reproduce Figure 3 from 1967 Hirsch Inertial Electrostatic Confinement
% Fusion paper.
%   - BUT use values of K+ and lambda+ from figure 5 to calculate 
%     ratio of real cathode radius to virtual anode radius
%   - Use figure 5 K+ & lambda+ values from the Erratum, not the original
%     paper
%
% RESULTS: 
%   1. K+=1.0 and lambda+=0.0 shows the real cathode at ~2.315 times
%      the virtual anode.
% 
%   2. K+=0.1 and lambda+=0.0 shows the real cathode at 3.98-4.26 times
%      the virtual anode, but potential never reaches close to 1????  Hmmm
%

% clean slate
clc; clear all; close all;

% radius from 1 to 4
max_radius = 3;
R = 1:0.01:max_radius;

% intial values at radius 1 (from the text):
%   v  = 0.0000001  - nudge slightly above zero to avoid infinity NaNs
%   v' = 0
hinits=[0.0000001,0];

K_plusses = [1,0.1];
lambda_plusses = [0,0];
legend_strs = cell(numel(lambda_plusses),1);

% several lambda_plus graphs on same graph
hold on

% axis labels and limits
title('Fig 3. Potential distributions');
ylabel('\phi(R)');
xlabel('R');
ax = gca;
ax.XTick = 0:1:max_radius;
ax.XLim = [0 max_radius];
ax.YTick = 0:0.2:1;
ax.YLim = [0 1.25];
ax.PlotBoxAspectRatio = [1,0.635,1];

grid on

% solving & plotting
label_idx = round((numel(R)-1)/3)+1;
for i=1:numel(lambda_plusses),
%for i=1:1,
    K_plus = K_plusses(i);
    lambda_plus = lambda_plusses(i);
    legend_strs{i} = sprintf('K_{+}=%0.4f, \\lambda_{+}=%0.3f',K_plus,lambda_plus);
    [R,v] = ode45(@hirsch_equ7,R,hinits,[],[K_plus lambda_plus]);
    plot(R,v(:,1)); % only plot the first column (second column is v')

    % create a silly text label
    alignment='left';
    label_str = sprintf('\\leftarrow %s',legend_strs{i});
    if logical(mod(i,2))
        alignment='right';
        label_str = sprintf('%s \\rightarrow',legend_strs{i});
    end
    %text(R(label_idx),v(label_idx),label_str,'HorizontalAlignment',alignment);
end

% add a legend
legend(legend_strs,'Location','southeast');

% no more plots for this figure
hold off

