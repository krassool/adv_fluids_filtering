function figure_format(fig_handle)

% Create axes
axes1 = axes('Parent',fig_handle,'FontSize',14,...
    'Position',[0.0551 0.0809 0.9185 0.873]);
box(axes1,'on');
hold(axes1,'on');

semilogx(X1,Y1); % DO PLOT
xlabel('y+','FontSize',14); % Create xlabel
ylabel('u-uinf/utau','FontSize',14); % Create ylabel

% Create title
title('Outer scaled velocity profile (deficit form)','FontSize',16);

