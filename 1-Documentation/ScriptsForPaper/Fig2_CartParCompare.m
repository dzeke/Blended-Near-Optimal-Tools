function Fig2_CartParCompare(fontsize,overplotpoint,axislabels)
% Generates a plot to compare cartesian and parallel coordinate plots sides by side
% This is Figure 2 in the citation below.
%
% EXAMPLE CALL
% Fig2_CartParCompae(14,1,{'Decision 1' 'Decision 2' sprintf('%s\n%s','Objective','Function')});
%
% INPUTS
%    fontsize = font size to plot labels and such
%
%    overplotpoint = boolean (1=yes) to over plot first point in a different
%            color. 0 to plot in same color
%
%    axislabels =  cell array of string for the x,y,z axis labels (default is X_1,
%             X_2, X_3). Note z is plotted as the first (left most) axis on
%             the parallel coordinate plot
%
%  CALLED FUNCTIONS
%   - \stats\stats\parallelcoords.m
%
% #####################
%   Programmed by David E. Rosenberg
%   2011 -- early use for Career proposal and conference proceedings papers
%   July 2014 -- added ability to label axes and position like for an
%       optimization problem with the objective function (z-axis) as the first
%       (left most) axis on the parallel coordinate plot
%      
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   Citation:
%   David E. Rosenberg (in review) "Near-optimal alternative generation,
%   visualization, and interaction for water resources decision making".
%   Water Resources Research. Submitted August 2014

%   Licensing:
%   This code is distributed AS-IS with no expressed or implied warranty regarding the claimed functionality. The entire code or parts 
%   may be used for any non-commercial purpose so long as the use is cited per the citation above. Use for any commercial purpose requires 
%   prior written permission from the author.
%
%   Bug Reports:
%   Bug reports are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, while much appreciated, there is no promise of when, if the bug will be corrected.

%%
% Define default axis labels if none provided
if nargin<3
    %Use default axis label names
    axislabels = {'X_1' 'X_2' 'X_3'};
end
    
%set up the data
X2 =  [3.5:-0.5:1]; 
X1 = [1:0.25:2.25];
X3 = X2 + X1; % A linear relationship

XFull = [X3' X1' X2'];

pcpOrder = [3 1 2];

Figure1 = figure;

%cartesian first
plot1 = subplot(1,2,1,'FontSize',fontsize-4);
if overplotpoint
   
   h = plot3(X1(2:end),X2(2:end),X3(2:end)); 
   hold on
   set(h,'linestyle','none','marker','o','color','b','markerfacecolor','b','markersize',10)
    
   h2 = plot3(X1(1),X2(1),X3(1)); 
   set(h2,'linestyle','none','marker','square','color',[1 0.526 1],'markerfacecolor',[1 0.526 1],'markersize',10);
   %add the text label
   text(X1(1),X2(1),X3(1),sprintf('(%.1f, %.1f, %.1f)',X1(1),X2(1),X3(1)),'fontsize',fontsize-2,'horizontalalignment','left','verticalalignment','bottom');
else
   h = plot3(X1,X2,X3); 
   hold on
   set(h,'linestyle','none','marker','o','color','b','markerfacecolor','b','markersize',10)

end

xlabel(axislabels{1},'fontsize',fontsize);
ylabel(axislabels{2},'fontsize',fontsize);
zlabel(axislabels{3},'fontsize',fontsize);
grid on
hold off

%now Parallel
plot2 = subplot(1,2,2,'FontSize',fontsize-4);
set(plot2,'Position',[.63 .3 .27 .55]);
hold on
h2 = parallelcoords(XFull,'parent',plot2);
set(h2,'linewidth',2,'color','b')
if overplotpoint
   h5 = parallelcoords(XFull(1,:),'parent',plot2,'linewidth',2,'color',[1 0.526 1]);
end

ymax = 5;
set(plot2,'YLim',[0 ymax],'Ytick',[0:1:5]);
ylabel('Value','fontsize',fontsize);

%Add text labels below each axis
% grab the y limits so we can position the axis labels
yLims = get(plot2,'YLim')

yPosAxisLabel = -0.04*(yLims(2)-yLims(1));

set(plot2, 'Xticklabel',[],'fontsize',fontsize-4)
for i=1:3
    text(i,yPosAxisLabel,axislabels{pcpOrder(i)},'fontsize',fontsize,'HorizontalAlignment', 'right','VerticalAlignment','middle','rotation',90);
end
%set(plot2,['x' 'TickLabel'], {'X_{1}' 'X_{2}' 'X_{3}'}); % label each tick with a string


%plot a fake x2 axis
h3 = plot([2 2],[0 ymax]);
set(h3,'linewidth',1,'color',[0 0 0])
h4 = plot([3 3],[0 ymax]);
set(h4,'linewidth',1,'color',[0 0 0])


