% set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendLocation','northwest','defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',14)
set(0, 'defaultTextFontSize',12);
set(0,'defaultLineLineWidth',1.5);


%{
%%
set(0,'DefaultFigureWindowStyle','normal')
i=1;
fig=figure(i)
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 12 12];
plot([0 1],[1 1])
print('5by3DimensionsFigure','-dpng','-r0')

%}