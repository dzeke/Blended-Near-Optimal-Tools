% ExampleTestPlots.m

% Test cases to illustrate syntax for using nearoptplotmo2.m
% and for debugging
%
% David E. Rosenberg
% November 2015
%
% -----------

ws = 1; %wait time between plots to allow screen to catchup

%% EXAMPLE USES
% A) Displays 11 evenly spaced traces across 5 axes (one objective function and 4 decision variables).
%    All settings at default values. 
%

nearoptplotmo2([0:.1:.9]',[[.1:.1:1]' [0:.1:.9]' [.1:.1:1]']);   
pause(ws)

% B) Displays 11 evenly spaced traces across 4 decision axes (no
%       objective funciton. Left ticks and scale same as right scale
nearoptplotmo2([],repmat([0:.1:1]',1,4));
pause(ws)

% C) Displays 15 traces across 5 axes with objective function values
%       determining the overall scale because they are largest. All other settings at default
%       values. Note the objective function scale obscures and "hides" the dynamic
%       ranges for the other axes.
nearoptplotmo2(150000*rand(15,1),[2+8*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'AxisScales','none');
pause(ws)

% C2) Improved plotting of Case C. Objective function plotted from 0 to
%       150,000 in units of 10^5 at left. Decision variables plotted in
%       units of 0, 2, 4, ... to 20 at right.
nearoptplotmo2(150000*rand(15,1),[2+8*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'AxisScales','custom', [1 2], [0 0 0 0 0; 150000 20 20 20 20]);
pause(ws)

% C3) Alternative plotting of Case C where minimum and maximum values for
% each axes are shown below and above each axis
nearoptplotmo2(150000*rand(15,1),[2+8*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'AxisScales','custom', [1 2], [0 2 0 0 0; 150000 10 20 5 5],'TickLabelPos','MinMaxEachAxis');
pause(ws)

% D) Labels the objective function and decision axes. Default AxisScales
%       setting gives similar scaling as Case C. Hides the panel of
%       controls. To make visible, uncheck Controls=>Hide all controls on
%       the menu.
nearoptplotmo2(150000*rand(15,1),[5*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'vObjLabels',{'Cost ($)'},'vXLabels',{'Area A' 'Area B' 'Machine 1' 'Machine B'},'AxisScales','auto','ShowControls',0);
pause(ws)

% E) Displays 30 traces across 9 axes. Assigns the first 20 rows to Group 1 plotted in light green and
%       last 10 rows to Group 2 plotted in blue.
%
nearoptplotmo2(150000*rand(30,2),[5*rand(30,1) 20*rand(30,1) 5*rand(30,6)],'vGroup',[ones(20,1);2*ones(10,1)],'AxisScales','auto','FontSize',18,'ShowObjsDiffColor',0);
pause(ws)

%       
% F) Like E but text labels used for groups. Note groups are ordered
%       alphanumerically so Group 1 still listed and plotted first in light green.
MyGroups = {'Group 2' 'Group 2' 'Group 2' 'Group 1' 'Group 1'}';
nearoptplotmo2(150000*rand(5,2),[5*rand(5,1) 20*rand(5,1) 5*rand(5,6)],'vGroup',MyGroups,'AxisScales','auto','FontSize',18,'ShowObjsDiffColor',0);
pause(ws)

% G) Show optimal solution (in thick black) and vertices (green) for a problem to maximize the objective cX subject to contraints that define the 5-dimensional simplex
%     of length 5. Once the plot loads, select the Interact Tab and click
%     the Generate button to see 150 sampled near-optimal alternatives.
    n = 5; r = 6; p = 150; % n = number of decision variables (dimensions), r=length of simplex, p=number of samples
    c = [5 3 1 0.5 0.5]; % c = objective function coefficients
    A = [-eye(n); ones(1,n)]; b = [zeros(n,1); r]; %Constraints of form Ax <= b define a simplex of length r
    mVert = [zeros(1,n); r*eye(n)]; %Vertices (corner points) of the simplex. Used for comparison on the plot.
   [Xopt,fopt] = linprog(-c,A,b);   %Solve for optimal solution, negative sign on c indicates maximization
   vObjsVert = mVert*c'; %Calculate objective function values for the vertices
   vGroups = ['Optimum'; repmat({'Vertices'},length(vObjsVert),1)];
   mGroupData = ['Vertices' {1} {1.5}; 'Optimum' {1} {2.5}];
   ProbForm.Aineq = [A;-c];
   ProbForm.bineq = [b;fopt];
   nearoptplotmo2([-fopt;vObjsVert],[Xopt';mVert],'Tolerance',0.85,'ProbForm',ProbForm,'cFunc',c,...
       'OptSolRow',1,'NearOptConstraint',7, 'ShowObjsDiffColor',0,'vGroup',vGroups,...
       'mGroupData',mGroupData,'GroupToHighlight','Optimum','GenerateType',3,'GenerateMethod',2,'NumSamples',p);  