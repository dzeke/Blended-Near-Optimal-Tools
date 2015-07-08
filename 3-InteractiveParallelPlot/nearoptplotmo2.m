function [hWindReturn] = nearoptplotmo2(mObjs, mDecisions, varargin)

% Makes a parallel coordinate plot and creates a graphical user interface (GUI) to allow a user to
% interact with the plot. The major functionalities are:
%
%  1) Option to link objective and decision spaces for an optimization problem and plot all on the same
%       parralel coordinate plot. (Objective Functions and their values [the
%       matrix mObjs] are plotted on the first 1..nO axes while
%       decision variable values in the matrix mDecisions are plotted on nD
%       subsequent axes to the right; pass an empty matrix for mObjs (mObjs=[]) to only show the mDecisions
%       matrix more akin to the standard parallelcoords function).
%
%       A row in each matrix mObjs and mDecisions represents a solution/alternative to the underlying optimization
%       problem (i.e., a "point" in the multivariate space).
%
%   2) Tools on the 'Display' tab (far right) to brush, pivot, highlight, and interactively manipulate the objective
%       function and decision variable data matricies as well as rows within them (i.e., both by columns/axes and rows/groups/color).
%
%   3) Tools on the 'Color Ramp' tab (middle) to select a particular axis
%       and ramp the color of lines of alternatives passing across the axis
%       according to value, e.g., from either dark to light color
%       (representing large to small values) or vice versa. This ramping
%       allows the viewer to simultaneously compare values on the checked
%       axis with values on each and every other axes.%       
%
%   4) Tools on the 'Interact' tab to generate new 
%       alterantives, e.g., for near-optimal analysis. There are three data sources available to use to generate
%       new alternatives: i) the data matrix itself (querying, i.e., a form of brushing), 2)
%       MATLAB (for linear programs using data specified in the objective
%       function vector(s) cFunc and a ProbForm structure that defines the inequality, equity, lower-, and upper-bound
%       contraints), and 3) executing a GAMS
%       file (note, you will need to install GAMS, see www.gams.com). The
%       methods can generate: a) one solution, b) the maximum extents
%       representing the minimum and maximum values given current selected
%       variables and their values, c) randomly generated alternatives, and d)
%       all alternatives via enumeration (for MIP problems)
%
%   5) Tools on the 'Plot Data' menu to save the current figure settings to the Matlab base workspace and use
%       to recreate the figure. Also to load additional data and print the current plot in
%       high-resolution PDF format.
%
%   6) Tools on the Control menu to show/hide various sets of controls on
%       the plot. These tools include all the controls, sliders, axes
%       checkboxes, alternatives to highlight, etc. 
%
%   7) Tools on the Update Formulation menu to streamline the process to add one or more new optimization model objective functions
%       or constraints and generate new alternatives from the updated
%       formulation (available only if the underlying problem is linear and
%       the existing objective function and constraint coefficients are
%       passed as optional parameters cFunc and ProbForm -- see #8 and
%       further description of optional inputs below)
%
%   8) Optional parameters/settings for each of the above features #2-7 to control how the parallel
%       coordinate plot is displayed and labeled and labeled upon run (input parameter varargin)
%
% OUTPUTS
%   hWindReturn - handle to plot generated. Will also return this value as
%   the parateter hWind if such a parameter is specified in varargin(see options below)

%
% INPUTS
%   mObjs - an m by nO matrix of objective function values (m solutions by
%      nO objectives). Leave as empty [] to only plot the mDecision matrix
% 
%   mDecisions - an m by nD matrix of decision variable values where each
%      row correspondings to the same row in mObjs (m solutions by nD decisions/variates)
%
%   varargin - an optional variable argument list comprising label-value
%      pairs {paramter1, value1, parameter2, value2, ...}
%      so that the user can specify optional parameter settings. Each pair starts with the parameter name as a text string
%      and is followed by the parameter value. Parameter names are case insensitive.
%      If a parameter is omitted or in incorrect format, the default value
%      is used. e.g., varargin = {'hWind',1,'AllowableDeviation',0.2} would set the
%      parameters hWind and AllowableDeviation, to the values, respectively,
%      of 1 and 0.2. Below is a listing of the parameters followed
%      by their meanings and default values (listed in parenthesis).
%
%     hWind = handle for figure window in which to plot. (Default: 0 = plots in a
%         new figure)
%
%     Tolerance = Near optimal tolerable deviation parameter. Fraction of the optimal objective function value and used to generate
%         new alternatives for optimization problems. Enter Tolerance > 1 for minimization objectives and Tolerance between 0 and 1 for maximization objectives.
%         Populates a text box on the Interact tab for generating additional data sets. For multiobjective problems this is a
%         (1 x nO) vector of values (Default value: 1 by nO vector of 1s)
%
%     AllowableDeviation = the deviation from the slider set value allowed when generating new alternatives
%         or highlighting alternatives on the plot. This value is specified in
%         units of fraction of a tick and can
%         also be dynamically changed by entering a value in the textbox on
%         the Interact Tab. The value allows the user to define a family of alternatives all having
%         values for the specified decision variables within the (Slider Set Value) +/- (Allowable Deviation)*(Tick Spacing)
%         (Default is 0 [no deviation]).
%         E.g., if the ticks for the axis are in units of 0, 5000, 10000,
%         ..., entering a AllowableDeviation of 0.1 will select all
%         alternatives within +/- 500 of the set value.
%
%     ErorrResid = the residual allowed when determining whether a
%           constraint is satisfied during Stratified Gibbs sampling.
%           (Default value: 0)
%
%     GenerateType = index that takes the values 1 (single solution), 2 (maximum extents), 3 (random sample), or 4 (enumerate for MIPs) to indicate
%           the type of new alternative(s) genreated when clicking the
%           Generate button on the Interact tab. Setting also available for user input as a drop-down list on the Interact tab. (Default: 1 [single
%           solution])
%
%     GenerateMethod = index that takes the values 1 (existing data set), 2 (Matlab linprog), or 3 (GAMS) to indicate
%           the computational method used to generate new alternative(s) genreated when clicking the
%           Generate button on the Interact tab. Setting also available for user input as a drop-down list on the Interact tab. (Default: 1 [existing data])
%
%     NumSamples = number of random samples to add when randomly sampling
%          additional poionts within the polytope defined by AMat and Brhs
%          (Default value = 0)
%
%     ObjSamplesPercent = the percent of samples to concentrate and
%          stratify along objective function axes/strata when stratify
%          sampling. Remaining samples are drawn from the decision variable
%          axes. (Default value= 15%)
%
%     sGamsFile = string with path and filename of GASMS file to run
%           that can generate new solutions to bring into the analysis (Default value: '') 


%     Precision = number of decimal points to which to truncate/round data
%           provided in mObjs and mDecisions to ensure uniformity of
%           plottin/rendering tools (default values: 2). Use negative values to
%           round to values larger than zero.
%
%     vFixed = vector with nO+nD elements that indicates
%           the axes value is fixed at the value vFixedVals [1=fixed,
%           0=no] (see explanation of vFixedVals below). This parameter can also be dynamically changed by
%           checking/unchecking one/more check boxes just above the
%           objective and decision axes labels. (default value:
%           all zeros [unchecked])

%     vFixedVals = vector of nO+nD elements whose values specify the
%           value at which to fix/set axes. The sliders, when visible, will also position at
%           these values. The user can also dynamically set these values by adjusting the sliders and/or entering
%           numeric values in the textboxes for the axes on the Sliders Tab or by clicking the Set All buttons on the same tab. (Default values: all zeros)

%     vStep = vector of nO+nD elements which specify the increments by which the sliders on the axes will move when the up/down arrows are clicked. 
%          In units of the decision variable value increments. (Default values: 1/20th the difference
%          between the minimum and maximum values in the dataset for the
%          axis).

%     FontSize = fontsize to display axis labels and all associated
%            deravative labels, controls, etc. (Default value: 16)
% 
%     mActCat = an nD by nF matrix whose rows represent the decision variables and nF columns are grouping fields to categorically classify
%           the decision variables; used to generate an optional table of decision variable labels to position below the decision variable axis labels, to
%           color-code groups of those labels/axes, and reorder axes (by groups) on the Display Tab (default value - [] or all ones)

%     vGroup  = m-length column vector whose values specify the group to
%           which each row in mObjs and mDecisions belongs and should be
%           colored on the plot. Can either be a numeric or text string. The unique values of vGroup are used to generate
%           the group labels and checkboxes on the Display tab where the user can dynamically specify which groups to show on the plot
%           (Default value: all ones [1 group])
%
%     vGroupOrder = an nU length column vector whose entires (unique values in vGroup) specify the
%           order in which to plot groups. The first group listed is plotted as the lower-most
%           layer and will appear first in the group list on the Display
%           tab (Default value: alpha-numeric sort)
%
%
%     mGroupData = an nU by 2 cell matrix whose entries describe attributes of each unique group in vGroup. 
%         Column 1 - vGroupName - The name of the group. The first group listed is plotted as the lower-most
%                   layer and will appear first in the group list on the
%                   Display Tab.
%         Column 2 - ShowGroup - binary variable (1=yes) to show the group
%                   on the plot and have the checkbox next to the group
%                   name checked.
%         Column 3 = Line Thickness. Thickness to plot lines associated
%                       with the group.
%       (Default values -- alphanumeric sort, all groups checked, and line thickness to 1)
%    
%     GroupToHighlight = value in vGroup representing a group to highlight on
%           the plot with thick lines of the color mHighlighColor (see below). E.g., use this feature to 
%           highlight the optimal solution (Default: 0 [no group])

%     CurrentRecord = row in mObjs and mDecisions to highlight on the plot.
%           (Default: 1 or the record number of a record in
%           GroupToHighlight if a GroupToHighlight is specified)%     
%
%     ShowCurrRecord = boolean takes the value of 1 to indicate to show/highlight the
%           current alternative selected on the Interact Tab. The feature
%           is also settable through the Controls=>Alternatives=>Show
%           current menu option (check to show) (Default: 0 [don't show])
%
%     ShowFilteredAlts = boolean takes the value of 1 to indicate to show/highlight 
%           filtered alternatives that match the set value critieria on checked axes. The feature
%           is also settable through the Controls=>Alternatives=>Show
%           filtered menu option (check to show) (Default: 0 [don't show])
%
%     mColors = an nU by 2 by d by 3 matrix indicating the color ramps to use
%           for the lines on the plot corresponding to each group. The first column indicates the
%           group. the second column is a binary value to indicate
%           either the decision axes (1) or objective axes (2). The third
%           column is the color depth (or hue) to use to plot the base
%           lines, highlighted subspaces, or moused-over solutions. The
%           fourth column contains the R-G-B values themselves. The default
%           color settings are OSU's Green (decision) and Magenta
%           (decision) hues from the 16 Step ramp
%           for the first group with darker hues representing sub-spaces
%           and moused over solutions. Then browns, blues, purples, reds,
%           oranges, and yellows (from the Categorical and
%           Stepped-Sequential OSU ramps) for subsequent groups.
%
%     mColorsAxisLabels = an nC*nF by 3 matrix indicating the colors to
%           use for axis labels and the group table lables. nG = number of decision variable categories
%           nF is the number of fields (grouping levels) within each
%           category. The axis label itself will be plotted in the color on
%           the nF, 2*nF, 3*nF, etc. rows. The final three columns
%           represent RGB values (Default is OSU Stepped Sequential scheme)
%
%     mColorYScales = a 2 x 3 matrix indicating the colors to use for the
%           ticks, Y-axis scales, and labels on the far right (row 1) and left (row 2) of the plot.
%           (Default values are color specified in mColor(1,1:2,3,:))
%
%     mHighlightColor = a 1 by 3 matrix indicating the colors to use for
%           highlighted lines on the plot (Group values =
%           GroupToHightlight). (Default value is Black [0 0 0])
%    
%     ShowControls = boolean takes the value of 1 to show the plot
%           controls as a panel of controls on the right of the main plot. Can
%           also be set from the Controls=>Hide all controls menu option.
%           (Default Value: 1 [show controls])
%
%     ShowGroupLabels = boolean takes the value of 1 to show a table below
%           the parallel coordinate axes labels that groups axes. 0 To high. Can
%           also be set from the Controls=>Show group labels menu option.
%           (Default Value: 0 [hide group labels]).
%
%     ShowLegend =  Controls whether to show the legend of traces. Three options:
%           OnAll - show all the time regardless of the state of
%                   ShowControls
%           OnHideControls - only show when ShowControls is set 0 (hide
%               controls)
%           Off - off all the time regardless of state.
%           Default value: OnHideControls. Note the legend will not show
%           when the Pareto Inset is shown.

%     ShowObjsDiffColor = boolean with a value of 1 to show the
%           objective function axes values, if they exist, in a different color
%           [specified in mColor(:,2,:,:)]. 0 to show the objective function axes
%           values in the same colors the mDecisions [i.e., mColor(:,1,:,:)].
%           Also available for user selection from the Controls menu.
%           (Default value: 0 [different color]).
%
%     ShowInsetPlot = boolean takes a value of 1 to show an inset plot in
%           cartesian coordinates for the 2 objectives. 0 hides. Ignores
%           if the problem does not have two objectives. (Default value: 0
%           [hide]).
%  
%     AxesScales = Parameter that determines how to transform columns of values in the mObjs and mDecisions matrices
%             so they can be plotted side-by-side on the same plot to see the dynamic ranges of values on each axis. The transformation
%             also updates the tick measurements reported on the left and right scales and possibly the axis labels. The AxisScales settings are:

%             1) none (DEFAULT) - no scaling. Plot mObjs and mDecisions in their native (original) units. The scale for each axis is the global scale
%                   generated by the maximum and minimum values in the dataset
%
%             2) standarize - akin to the Standarize option on ParrallelCoorids where
%                   zero is the mean for data in the column and 1/-1 are one standard
%                   deviation above/below the mean
%
%             3) normalize - Transforms the data so all values are plotted between 0 (repreenting the minimum data value for
%                   the column) and 1 (representing the maximum value).
%
%             4) auto  - the program auto-calculates the transformations to display the ranges of values on each axis in convienent
%                   tick mark increments of generally 1, 2, 3, 5, or 10s.
%                   The scaling will simulataneously update and note the
%                   magnitude of each axis as part of the label for the
%                   decision or objective axis.
%
%             5) custom - generates a custom transformation for each axis specified by additional lower and upper bound limits on each axis the user provides
%                   as well as which objective and decision axes to use as the base for showing tick marks. The additional parameters are BaseAxes and mLims and can be specified
%                   in either the form {'AxisScales','custom',BaseAxes,mLims} or {'AxisScales','custom','BaseAxis',BaseAxes,'mLims',mLims}. 
%
%                   BaseAxes -- a 2 element vector whose values represent
%                       the index of the objective and decision axes to use as
%                       the base scales for showing ticks on the left and right. i.e. a setting of [1 3] means the
%                       first obective axis will be used as the scale and ticks at the plot left to
%                       show all objective function axes (subsequent obj
%                       axes values will be transformed to that scale). 
%                       Similarly, the 3rd decision variable axis will be
%                       used as the scale and ticks to show all decision variable
%                       axes at the far right (all other decision variable
%                       axes will be transformed to that scale)
%
%                   mLims is a 2 by (nO+nD) matrix of user-supplied lower (row 1) and upper (row 2) limits
%                         for each axis. Using the user-supplied limits,
%                         the program transforms the data for the axis into
%                         values comensorate with Base Axis. The transformation is appended to the bottom of the decision label. Does 
%                         separate transformations for objective axes and
%                         decision axes. E.g., a setting of [0 0 0 0 0;
%                         1000 20000 20000 200 2000000] for a dataset with
%                         1 objective function, 4 decision variables, and
%                         the BasesAxes setting [1 3] will scale the
%                         objective function scale between 0 and 1,000 and
%                         the decision axis scale from 0 to 200. The first
%                         and second decision axes will be transformed so
%                         that the values 0 to 20,000 plot at 0 to 200, and
%                         the fourth axis so values 0 to 2,000,000 plot
%                         between 0 and 200. The labels for axes 1, 2, and
%                         4 will include the multipliers 100, 100, and
%                         10,000 to indicate the values shown on the plot
%                         need to be mutliplied by those factors to get the
%                         data value.

%`           Once a plot is generated, hover the mouse over an axis to see
%               the actual data value.
%
%     NumTicks = specifies the number of ticks to include on the left and
%               scales. A scalar value indicates the same number of ticks
%               on the left and right scales. A 2 element vector, e.g, 
%               [8, 4] specifies 8 ticks on the left scale and 4 ticks on
%               the right scale. NOTE: the right number of ticks is ignored for
%               AxisScales types of Normalize, Standardize, and None (these plots require
%               the same number of ticks on the left and right)(default value 6 -- 6 ticks on the left
%               and right scales)
%
%     TickMag = a scalar that specifies the mangitude (in base 10 units) by
%              at which to truncate zeros in the tick labels and instead
%              show them as part of the axis label. e.g., a TickMag of 4 would cause
%              a tick value of 1,000,000 to be shown as 100 and the y axis
%              label to read "Value (10^4)". (Default value = 5)
%     
%     ProbForm = a Matlab structure for the linear programming formulation that defines the
%           inequality, equity, lower-bound, and upper-bound constraints for the problem
%           as the fields
%               .Aineq = lCon x n matrix of inequality constraint coefficients
%               .bineq = mlCon x 1 vector of right-hand side constraint coefficients for inequality constraints
%               .Aeq  = q x n matrix of equality constraint coefficients          
%               .beq  = q x 1 vector of right-hand side constraint
%                       coefficients for equality constraints
%               .lb   = n x 1 vector of lower-bound limits for decision
%                           variables
%               .ub   = n x 1 vector of upper-bound limits for decision
%                       variables
%            (i.e., .Aineq x X <= .bineq ; .Aeq x X = .beq; X >= .lb; X <= .ub)
%            Other fields of ProbForm used as is, but do not pass objective
%            function as f. Rather use cFunc (see below).
%
%               .solver = solver to use solve the model (Default: linprog)
%               .options = options to pass to the solver (Default:
%               structure with interior-point algorithm, maxiter 35000,
%               display off
%
%           (Default value [] for structure or [] for fields)
%
%     AMat = OLD WAY OF HANDLING CONSTRAINTS lCon by nD matrix of coeffiecients for the lCon constraints
%           in the underlying linear optimization problem associated with the
%           data provided. This parameter, along with Brhs and possibly cFunc, is used when clicking the 'Generate' button on the
%           Interact tab to sample new alternatives. Constraints need to be defined in
%           form AMat * X <= Brhs (default value [])
%
%     Brhs = OLD EAY OF HANDLING CONSTRAINTS column vector of lCon elements whose values are the right hand side
%           of the constraints that define the system of equations from which to generate additional alternatives.
%           i.e., AMat * X <= Brhs. (Default value [])
%
%     cFunc = nO by nD matrix of coeffiecents where each row respresents
%           the coefficients for an objective function of the
%           optimization problem. Used with AMat and Brhs to calculate the
%           resulting objective function values for newly generated alternatives.
%           The underlying optimization problem for the first
%           objective (row) in cFunct is defined as:
%              Min Z1 = cFunc(1,:)*X
%               subject to: constraints in ProbForm
%           (Default value [])
%
%`    NearOptConstraint = Index that references a row in ProbForm.Aineq that is the
%           near-optimal tolerable deviation constraint. For a multi-objective problem, this
%           will be a (1 x nO) vector of indexes. Together with cFunc and OptSolRow
%           (below) this parameter is used to update the near-optimal constraint
%            when generating alternatives. (Default value []).
%
%     OptSolRow = Index that references a row in mObjs and mDecisions that
%           is the optimal solution. For a multi-objective problem, this
%           will be a (1 x nO) vector of indexes. Together with cFunc and
%           NearOptConstraint (above) this parameter is used to update the near-optimal constraint
%            when generating alternatives.  (Default value []).
%
%     StartTab = integer or cell that specifies the tab index to show as
%           visible on start up. User can click tab buttons to dynamically change. 
%           Used to "remember" the most recent selected tab 
%           when re-plotting. (Default Value: 3 [Display])
%
%     HideSliders = boolean with value of one that means the Hide Sliders
%           on the interact tab checkbox is checked and sliders are not shown
%           when the plot opens. 0 = show sliders. User can click the checkbox to dynammically change.
%           Used to "remember" the most recent setting
%           when re-plotting (Default: 1 [sliders hidden])
%
%     PlotPosition = position of plot within the figure in standard
%           normalized units [left bottom width height]. The width leaves room
%           at the right for the control tabs.
%           (Default of [0.100 0.47 .615 0.6182 - yBottom + 0.3068])
%
%     YAxisMargins = 1 x 2 vector that specifies the margin between the left 
%           YAxis ticks and labels and the left-most parallel axis (1st value) and
%           and the right-most parallel axes and right YAxis ticks and labels.
%           Expressed as a fraction of the parallel axis spacing. Default
%           value: no spacing on left and 1/2 the axes spacing on the right [0 0.5].
%           A single value gives the same margin on the left and right.
%
%     PanelWidth = width of panel of controls at the far right in
%           pixels. This width will overide the width setting in
%           PlotPostion(3). Use to ensure the controls are visible.
%           (Default value [285])
%
%     vObjLabels = row vector of nO text labels for the objective function axes.
%          (Default values are {'Obj 1' 'Obj 2' ...})
%
%     vXLabels = 1 x nD (row) vector of nD text labels for decision variable axes. (Default values are {'Decision 1' 'Decision 2' ...})
%
%     vXLabelsShort = 1 x nD (row) vector short text labels for decision
%           variables. Used should the user decide to recalculate and load additional results
%           from GAMS (default value is same as vXLabels)
%
%     yAxisLabels = 1 x 2 cell array of string texts that indicate the labels for the {left-axis
%           right-axis}. If you specify units in parenthesis, e.g., 'Left
%           Label (Ac-Ft)', the program may add a scalar in front if an
%           AxisScales option is selected. (Default value = {'Objective
%           Functions' 'Decision Variables'})
%
%
%
% EXAMPLE USES
%
% A) nearoptplotmo2([0:.1:.9]',[[.1:.1:1]' [0:.1:.9]' [.1:.1:1]']);
%       Displays 11 evenly spaced traces across 5 axes (one objective function and 4 decision variables). All settings at default values. 
%
% B) nearoptplotmo2([],repmat([0:.1:1]',1,4));
%       Displays 11 evenly spaced traces across 4 decision axes (no
%       objective funciton. Left ticks and scale same as right scale
%
% C) nearoptplotmo2(150000*rand(15,1),[5*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'AxisScales','none');
%       Displays 15 traces across 5 axes with objective function values
%       determining the overall scale because they are largest. All other settings at default
%       values. Note the objective function scale obscures and "hides" the dynamic
%       ranges for the other axes.
%
% C2) nearoptplotmo2(150000*rand(15,1),[5*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'AxisScales','custom', [1 2], [0 0 0 0 0; 150000 20 20 20 20]);
%       Improved plotting of Case C. Objective function plotted from 0 to
%       150,000 in units of 10^5 at left. Decision variables plotted in
%       units of 0, 2, 4, ... to 20 at right.
%
% D) nearoptplotmo2(150000*rand(15,1),[5*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'vObjLabels',{'Cost ($)'},'vXLabels',{'Area A' 'Area B' 'Machine 1' 'Machine B'},'AxisScales','auto','ShowControls',0);
%       Labels the objective function and decision axes. Default AxisScales
%       setting gives similar scaling as Case C. Hides the panel of
%       controls. To make visible, uncheck Controls=>Hide all controls on
%       the menu.
%
% E) nearoptplotmo2(150000*rand(30,2),[5*rand(30,1) 20*rand(30,1) 5*rand(30,6)],'vGroup',[ones(20,1);2*ones(10,1)],'AxisScales','auto','FontSize',18,'ShowObjsDiffColor',0);
%       Displays 30 traces across 9 axes. Assigns the first 20 rows to Group 1 plotted in light green and
%       last 10 rows to Group 2 plotted in blue.
%
% F) MyGroups = {'Group 2' 'Group 2' 'Group 2' 'Group 1' 'Group 1'}';
%    nearoptplotmo2(150000*rand(5,2),[5*rand(5,1) 20*rand(5,1) 5*rand(5,6)],'vGroup',MyGroups,'AxisScales','auto','FontSize',18,'ShowObjsDiffColor',0);
%       Like E but text labels used for groups. Note groups are ordered
%       alphanumerically so Group 1 still listed and plotted first in light green.
%
% G) n = 5; r = 6; p = 150;
%    c = [5 3 1 0.5 0.5];
%    A = [-eye(n); ones(1,n)]; b = [zeros(n,1); r]; %Constraints of form Ax <= b define a simplex of length 6                                                 
%    mVert = [zeros(1,n); r*eye(n)]; %Vertices (corner points) of the simplex. Used for comparison on the plot.
%   [Xopt,fopt] = linprog(-c,A,b);   %Solve for optimal solution, negative sign on c indicates maximization
%   vObjsVert = mVert*c'; %Calculate objective function values for the vertices
%   vGroups = ['Optimum'; repmat({'Vertices'},length(vObjsVert),1)];
%   mGroupData = ['Vertices' {1} {1.5}; 'Optimum' {1} {2.5}];
%   ProbForm.Aineq = [A;-c];
%   ProbForm.bineq = [b;fopt];
%   nearoptplotmo2([-fopt;vObjsVert],[Xopt';mVert],'Tolerance',0.85,'ProbForm',ProbForm,'cFunc',c,...
%       'OptSolRow',1,'NearOptConstraint',7, 'ShowObjsDiffColor',0,'vGroup',vGroups,...
%       'mGroupData',mGroupData,'GroupToHighlight','Optimum','GenerateType',3,'GenerateMethod',2,'NumSamples',p);
%
%     Show optimal solution (in thick black) and vertices (green) for a problem to maximize the objective cX subject to contraints that define the 5-dimensional simplex
%     of length 5. Once the plot loads, select the Interact Tab and click
%     the Generate button to see 150 sampled near-optimal alternatives.
%   
%% #####################
%   Programmed by David E. Rosenberg
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%   - July 2012 First version to include multi-objective problems
%   - Updated May 2013 to include passing optional input arguements as a variable list and
%      using the sliders and axes checkboxes to view subsets of solutions
%   - Updated June-July 2014 to include alternative generation and interaction functions, reorganize GUI, allow data save, 
%      and improve documentation.
%   - Updated Sept 2014 to add Update Formulation menu -- further improve the process to streamline adding
%      additional objectives and constraints
%   - Updated Feb 2015 to define model formulation with ModForm Matlab
%      structure to allow inequality, equality, lower-bound, and
%      upper-bound constraints
%
%   Citation:
%   David E. Rosenberg (2015) "Blended near-optimal alternative generation, 
%   visualization, and interaction for water resources decision making".
%   Water Resources Research. doi:10.1002/2013WR014667.
%   http://onlinelibrary.wiley.com/doi/10.1002/2013WR014667/full
%
%   Licensing:
%   This code is distributed AS-IS with no expressed or implied warranty regarding functionality. The entire code or parts 
%   may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires 
%   prior written permission from the author.
%
%   Bug Reports and Feedback:
%   This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.
%
%% #######################

% Read in stuff about the main input parameters
%  
    [m,nD]=size(mDecisions);
    
    if isempty(mObjs)
        nO = 0;
        mResults = mDecisions;
    else
        [mO,nO] = size(mObjs);

        if (mO~=m)
            error('nearoptplot2:bad dataset',sprintf('%.f rows of mObj do not match %.f rows of mDecisions',mO,m));
        end
        %Aggregate the Objective and Decision matricies
        mResults = [mObjs mDecisions];
    end   
      
    n = nO+nD;
    
    %% Set default field values
    hWind = 0;
    Tolerance = ones(1,nO);
    AllowableDeviation = 0; ErrorResid = 0;
    NumSamples=0;
    ObjSamplesPercent = 15; %Percent of samples drawn from objective function axes/strata
    Precision=2;
    vFixed = zeros(n,1);
    vFixedVals = vFixed;
    vStep = NaN; %dummy setting to later use to set to default value of if no value passed; (max(mResults)-min(mResults))/20;
    
    vObjLabels = cell(1,nO);
    for i=1:nO
        vObjLabels{i} = ['Objective ', num2str(i)];
    end
    vXLabels = cell(1,nD);
    for i=1:nD
         vXLabels{1,i} = ['Decision ', num2str(i)];
    end
    
    vXLabelsShort = vXLabels;
    yAxisLabels = {'Objective Functions' 'Decision Variables'};

    FontSize = 16;
    mActCat = repmat({'1'},nD,1); lShowGroupLabelsEnable = 'off';
    AxisScales='none';
    BaseAxis=[min(nO,1) 1];  %put a zero in the first element if there are no objective axes
    NumTicks = [6 6];
    TickMag = 5;
    mLims = zeros(2,n);
    
    iOpt = '0';
    lCurrRecord = 1; %First record
    lShowCurrRecord = 0; %Do not show
    lShowFiltered = 0; %Do not show
    
    vGroup = ones(m,1);
    vGroupOrder = [];
    mGroupData = {};
    mColors = []; mColorsAxisLabels = []; lShowObjsDiffColor = 1; 
    mColorYScales = [];
    mHighlightColor=[0 0 0];
    
    vMaxValues = max(mDecisions);
    ProbForm = []; ProbFormTemp=[];
    AMat = []; AMatTemp=[];
    Brhs = []; BrhsTemp=[];   
    cFunc = []; cFuncTemp=[];
    NearOptConstraint = []; NearOptConstraintTemp = [];
    OptSolRow = []; OptSolRowTemp = [];
    blValidLinProgramData = 0; %Assume for starters there is not valid linear programming data in AMat, Brhs, cFunc
    
    HideSliders = 1; HideCheckboxes = 0;
    GenType = 1;
    GenMethod = 1;
    sGamsFile = '';
    
    
    yBottom = 0.47;
    PlotPosition = [0.100 yBottom 0.615 0.6182 - yBottom + 0.3068];
    YAxisMargins = [0 0.5];
    PanelWidth = 60; %Characters; PanelWidth = 285; %pixel
    
    cLegendOptions = {'Off','OnHideControls','OnAll'};
    ShowControls = 1; lShowGroupLabels = 0; ShowInsetPlot = 0; ShowLegend = cLegendOptions{2};
    
    ButtonText = {'Interact' 'Color Ramp' 'Display'};
    StartTab = 3;
    lStartTab = StartTab; %Integer value
    cStartTab = ButtonText(lStartTab); %Cell value
    
    %Fields to ignore
    IgnoreFields = {'mConvert' 'AllowableDeviationPU' 'nO', 'SubSpaceGroup' 'SubSpaceError' 'SubSpaceErrorPU'};
    
    %% Read in the optional parameters from the variable argument list
    % to overide default values
    if nargin >= 3
        count=0;
        varargin;
        PrevCnt=0;
        
        while (count<nargin-2),
            count = count+1;
            StCount = count;
            
            if (ischar(varargin{count}) && sum(strcmpi(varargin{count},IgnoreFields))>0)
                %Ignore the field
                count=count+1;
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'hWind'))
                hWind = varargin{count+1};
                count=count+1;
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'Tolerance'))            
                TempTolerance = varargin{count+1};
                
                if (nO==0) % ignore the input; continue with empty
                    
                elseif (size(TempTolerance,2) ~= nO) || (size(TempTolerance,1) ~= 1)
                    warning(['nearoptplotmo2: Tolerance must be a single row vector with the same number of columns as mObjs (', num2str(nO), '). Continuing with default of all 1s.'])
                else
                    Tolerance = TempTolerance;
                end
                
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'AllowableDeviation'))            
                AllowableDeviation = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'ErrorResid'))            
                ErrorResid = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'NumSamples'))            
                NumSamples = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'ObjSamplesPercent'))            
                ObjSamplesPercent = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'GenerateType'))            
                GenType = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'GenerateMethod'))            
                GenMethod = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'HideSliders'))            
                HideSliders = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'HideCheckboxes'))            
                HideCheckboxes = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'Precision'))            
                Precision = varargin{count+1};
                count=count+1;
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vFixed'))
                if (length(varargin{count+1}) == nO+nD)
                    vFixed = varargin{count+1};
                else
                    warning(['nearoptplotmo2: vFixed is of different size than ', num2str(nO+nD), '. Continuing with default.'])
                end
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vFixedVals'))
                if (length(varargin{count+1}) == nO+nD)
                    vFixedVals = varargin{count+1};
                else
                    warning(['nearoptplotmo2: vFixedVals is of different size than ', num2str(nO+nD), '. Continuing with default.'])
                end
                count=count+1;
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vStep'))
                if (length(varargin{count+1}) == n)
                    vStep = varargin{count+1};
                else
                    warning(['nearoptplotmo2: vStep is of different size than ', num2str(n), '. Continuing with default.'])
                end
                count=count+1;
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mActCat'))
                
                if (size(varargin{count+1},1) == nD)
                    mActCat = varargin{count+1};
                    lShowGroupLabelsEnable = 'on';
                else
                    warning(['nearoptplotmo2: mActCat has a different number of rows than the number of columns in mDecisions', num2str(nD), '. Ignoring input.'])
                end
                count=count+1;

            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'AxisScales'))
                
                AxisScales = lower(varargin{count+1});
                
                switch (AxisScales)
                    case 'none'
                        count=count+1;
                    case 'standardize'
                        %make the plot window narrower to accomate the
                        %extra tick labels
                        %PlotPosition(1) = .15;
                        %PlotPosition(3) = PlotPosition(3)-.1;
                        count=count+1;
                    case 'normalize'
                        count=count+1;
                    case 'custom'
                        if (nargin>=count+2) && (isnumeric(varargin{count+2}))
                            BaseAxis = varargin{count+2};
                            count=count+2;
                        else
                            count=count+1;
                        end
                        if (nargin>=count+1) && (isnumeric(varargin{count+1}))
                            mLims = varargin{count+1};
                            count=count+1;
                        end

                    otherwise
                        if ~strcmpi(AxisScales,'auto')
                            warning(['nearoptplotmo2: AxisScale value of ', AxisScales, ' not recognized. Continuing with default setting of auto.'])
                        end

                        AxisScales='auto';
                                                   
                        count=count+1;
                end
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'BaseAxis'))                
                BaseAxis = varargin{count+1};    
                count=count+1;

            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mLims'))                
                mLims = varargin{count+1};    
                count=count+1;                                
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vGroup'))
                varargin{count};
                varargin{count+1};
               
                if (size(varargin{count+1},1) == m)
                   vGroup = varargin{count+1};
                else
                    warning(['nearoptplotmo2: vGroup has a different number of rows than the data set', num2str(m), '. Continuing with default of one group.'])
                end
                count=count+1;
                
             elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vGroupOrder'))             
                vGroupOrder = varargin{count+1};
                count=count+1; 
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mGroupData'))             
                mGroupData = varargin{count+1};
                count=count+1; 
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'FontSize'))
                FontSize = abs(varargin{count+1});
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'NumTicks'))
                NumTicks = abs(varargin{count+1});
                
                if size(NumTicks,1) == 1
                    NumTicks = [NumTicks NumTicks];
                elseif (size(NumTicks,1) > 2)
                    NumTicks = NumTicks(1:2);
                    warning(['nearoptplotmo2: NumTicks has more than 2 elements. Continuing with first two element values.'])
                end

                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'TickMag'))
                TickMag = abs(varargin{count+1});
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'GroupToHighlight'))
                iOpt = varargin{count+1};                
                if isnumeric(iOpt)
                    iOpt = num2str(iOpt);
                end
                count=count+1;
 
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'CurrentRecord'))
                lCurrRecord = varargin{count+1};
                count=count+1;         
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'ShowCurrRecord'))
                lShowCurrRecord = varargin{count+1};
                count=count+1;   
                               
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'ShowFilteredAlts'))
                lShowFiltered = varargin{count+1};
                count=count+1;   
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mColors'))
                mColors = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mColorsAxisLabels'))
                mColorsAxisLabels = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mColorYScales'))
                mColorYScales = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mHighlightColor'))
                if (size(varargin{count+1},2) == 3)
                   mHighlightColor = varargin{count+1};
                else
                    warning(['nearoptplotmo2: mHighlightColor must only be a 1 by 3 matrix. Continuing with default highlight color setting.'])
                end
                
                count=count+1;
                
             elseif (ischar(varargin{count}) && strcmpi(varargin{count},'ShowObjsDiffColor'))
                lShowObjsDiffColor = varargin{count+1};
                count=count+1;              
 
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'ProbForm'))
                ProbFormTemp = varargin{count+1};
                if ~strcmpi(class(ProbFormTemp),'struct')
                    warning('nearoptplotmo2: ProbForm is not a structure. Defaulting to no problem formulation');
                    ProbFormTemp = [];
                end
                count=count+1;
                
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'AMat'))
                AMatTemp = varargin{count+1};
                if strcmpi(AMatTemp,'cell')
                    AMatTemp = cell2mat(AMatTemp);
                end
                count=count+1;

            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'Brhs'))
                BrhsTemp = varargin{count+1};
                if strcmpi(BrhsTemp,'cell')
                    BrhsTemp = cell2mat(BrhsTemp);
                end
                count=count+1;
            
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'cFunc'))
                cFuncTemp = varargin{count+1};
                if strcmpi(cFuncTemp,'cell')
                    cFunc = cell2mat(cFuncTemp);
                end
                
                [mCF nCF] = size(cFuncTemp);
                
                if ~isempty(cFuncTemp) && ((mCF~=nO) || (nCF~=nD))
                    warnmsg = sprintf('nearoptplotmo2: cFunc is improperly sized. Must be %d by %d. Default to empty [].',nO,nD);
                    cFuncTemp = [];
                    warning(warnmsg)
                end
                count=count+1;
                
           elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'NearOptConstraint'))
                NearOptConstraintTemp = varargin{count+1};
                if ~isempty(NearOptConstraintTemp) && ((size(NearOptConstraintTemp,2) ~= nO) || (size(NearOptConstraintTemp,1) ~= 1))
                    warnmsg = ['nearoptplotmo2: NearOptConstraint must be a single row vector with the same number of columns as as mObjs. Default to empty []'];
                    warning(warnmsg)
                else
                    NearOptConstraint = NearOptConstraintTemp;
                end
                count=count+1;            
                
                
           elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'OptSolRow'))
                OptSolRowTemp = varargin{count+1};             
                if any(OptSolRowTemp < 0) || any(OptSolRowTemp > m)
                    warnmsg = ['nearoptplotmo2: OptSolRow must be > 0 and less than ',num2str(m),'. Default to empty [].'];
                    warning(warnmsg)                   
                elseif ~isempty(OptSolRowTemp) && ((size(OptSolRowTemp,2) ~= nO) || (size(OptSolRowTemp,1) ~= 1))
                    warnmsg = ['nearoptplotmo2: OptSolRow must be a single row vector with the same number of columns as mObjs. Default to empty []'];
                    warning(warnmsg)
                else
                    OptSolRow = OptSolRowTemp;
                end
                
                count=count+1;                       
            
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'PlotPosition'))
                PlotPosition = varargin{count+1};
                
                if (length(varargin{count+1})==4)
                    PlotPosition = varargin{count+1};
                    yBottom = PlotPosition(2);
                else
                    warnmsg = ['nearoptplotmo2: PlotPosition has', num2str(length(PlotPosition)), ' but needs 4. Continuing with default position.'];
                    warning(warnmsg)
                end
                count=count+1;
                
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'YAxisMargins'))
                
                YAxisMarginsTemp = varargin{count+1};
                              
                if (isnumeric(varargin{count+1})) && (size(YAxisMarginsTemp,2) <= 2) && (YAxisMarginsTemp(1)>=0) && ((size(YAxisMarginsTemp,2) == 2) && (YAxisMarginsTemp(2)>=0))
                    if size(YAxisMarginsTemp,2)==1
                        YAxisMargins = [YAxisMarginsTemp YAxisMarginsTemp];
                    else
                        YAxisMargins = YAxisMarginsTemp;
                    end
                else
                    warnmsg = ['nearoptplotmo2: Error with YAxisMargins. Continuing with default value.'];
                    warning(warnmsg)
                end
                count=count+1;              
                
             elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'PanelWidth'))
                              
                if (isnumeric(varargin{count+1}))
                    PanelWidth = varargin{count+1};
                else
                    warnmsg = ['nearoptplotmo2: PanelWidth is not numeric. Continuing with default position.'];
                    warning(warnmsg)
                end
                count=count+1;

             elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'StartTab'))
                tmpStartTab = varargin{count+1};
                %handling depends on class
                if strcmpi(class(tmpStartTab),'cell') || strcmpi(class(tmpStartTab),'char')
                        cnt = [1:3];
                        ltmpStartTab = cnt(strcmpi(ButtonText,tmpStartTab));
                        blUseProvided = ~isempty(ltmpStartTab);
                        tmpStartTabError = tmpStartTab;
                        dfStartTab = cStartTab;
                        
                elseif strcmpi(class(tmpStartTab),'double')
                        ltmpStartTab = round(tmpStartTab);
                        tmpStartTabError = num2str(lStartTab);
                        blUseProvided =  (ltmpStartTab >= 1) && (ltmpStartTab <= 3);
                        dfStartTab = num2str(lStartTab);
                        
                else
                        blUseProvided = 0;
                end
                
                if blUseProvided
                    StartTab = tmpStartTab;
                    lStartTab = ltmpStartTab;
                    cStartTab = ButtonText(lStartTab);
                else
                    warnmsg = ['nearoptplotmo2: StartTab of ', tmpStartTabError, ' not recognized. Continuing with default of ', dfStartTab,'.'];
                    warning(warnmsg)                            
                end                   
                
                count=count+1;
                              
             elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'ShowControls'))
                ShowControls = varargin{count+1};
                %convert to numeric
                switch class(ShowControls)
                    case 'cell'
                        ShowControls = cell2mat(ShowControls);
                    case 'char'
                        ShowControls = str2num(ShowControls);
                end
                %reset PlotPosition if not showing controls
                if ShowControls == 0
                    PlotPosition = [0.100 yBottom 0.8 1-yBottom-0.1];
                end
                
                count=count+1;       
                
             elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'ShowLegend'))                
                ShowLegendTemp = varargin{count+1};
                if ismember(lower(ShowLegendTemp),lower(cLegendOptions))
                    ShowLegend=ShowLegendTemp;
                else
                    warning(['nearoptplotmo2: ShowLegend setting ', ShowLegendTemp, ' not recognized. Ignoring. Proceeding with default ',ShowLegend])
                end
                count=count+1;                      
                
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'ShowGroupLabels'))
                lShowGroupLabels = varargin{count+1};
                %convert to numeric
                switch class(lShowGroupLabels)
                    case 'cell'
                        lShowGroupLabels = cell2mat(lShowGroupLabels);
                    case 'char'
                        lShowGroupLabels = str2num(lShowGroupLabels);
                end
                
                count=count+1;    
                
                
            elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'ShowInsetPlot'))
                ShowInsetPlotTemp = varargin{count+1};
                
                if (nO==2)
                    ShowInsetPlot=ShowInsetPlotTemp;
                elseif (ShowInsetPlotTemp > 0)
                    warning(['nearoptplotmo2: Ignoring ShowInsetPlot input. ', num2str(nO), ' objectives provided. Must have 2 objectives to show Inset Plot. Continuing with default.'])
                end
                
                count=count+1;  
              
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vObjLabels'))
              %Objective Function Labels
                if (length(varargin{count+1}) == nO)
                    vObjLabels = varargin{count+1};
                else
                    warning(['nearoptplotmo2: vObjLabels is of different size than ', num2str(nO), '. Continuing with default.'])
                end
                
                count=count+1;
                
              elseif (ischar(varargin{count}) && strcmpi(varargin{count},'yAxisLabels'))
                % Y Axis Labels
                    if (length(varargin{count+1}) == 2)  
                        yAxisLabels = varargin{count+1};
                    elseif (nO==0) && (length(varargin{count+1}) == 1) 
                        yAxisLabels = {{} varargin{count+1}};
                    else
                        warning(['nearoptplotmo2: yAxisLabels is of different size than 2. Continuing with default.'])
                    end

                    count = count+1;

              elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vXLabels'))
                % Decision axis labels
                    if (length(varargin{count+1}) == nD)
                        vXLabels = varargin{count+1};
                    else
                        warning(['nearoptplotmo2: options.vXLabels is of different size than ', num2str(nD), '. Continuing with default.'])
                    end

                    count = count+1; 

               elseif (ischar(varargin{count}) && strcmpi(varargin{count},'vXLabelsShort'))
                % Decision axis labels
                    if (length(varargin{count+1}) == nD)
                        vXLabelsShort = varargin{count+1};
                    else
                        warning(['nearoptplotmo2: vXLabelsShort is of different size than ', num2str(nD), '. Continuing with default.'])
                    end

                    count = count+1;
                    
                elseif (ischar(varargin{count}) && strcmpi(varargin{count},'TickSize'))
                % TickSize, ignore this input
                    count = count+1;
                    
                    
                elseif (ischar(varargin{count}) && strcmpi(varargin{count},'sGamsFile'))
                % Gams File name
                                   
                    if (exist(varargin{count+1},'file'))
                        sGamsFile = varargin{count+1};
                    elseif ~strcmpi(varargin{count+1},'')
                        warning(['nearoptplotmo2: options.sGamsFile (%s) does not exist ', varargin{count+1}, '. Continuing with default.'])
                    end

                    count = count+1; 
 
              else
                  %do not recognize field
                  if ischar(varargin{count})
                      msg = ['Do not recognize parameter ',varargin{count}];
                  else
                      if PrevCnt==0
                          PrevField = '1st Parameter';
                      else
                          PrevField = ['Parameter after ',varargin{PrevCnt}];
                      end
                      
                      msg = [PrevField, ' is missing a parameter identifier'];   
                  end
                  error(['nearoptplotmo2: ', msg])
                    
            end
              PrevCnt = StCount;
        end
        
        %% Additional Error checking on some of the optional parameters
        
        %error check the BaseAxis and mLims parameters
        if (length(BaseAxis) ~= 2) || ((BaseAxis(1)>nO) && (nO>0)) || (BaseAxis(2)>nD)
            warning(['nearoptplotmo2: AxisScale BaseAxis does not have 2 elements or values are larger than the number of axes. Continuing with default values of 1st objective and decision axes.'])
            BaseAxis = ones(1,2);
        end
        if (size(mLims,2) ~= nO+nD) || (size(mLims,1)<2)
            warning(['nearoptplotmo2: AxisScale mLims is improperly sized. Continuing with default AxisScales setting of auto.'])
            AxisScales='auto';
        end
        count=count+3;
        
        %Further error checking on the ProbForm structure
        if ~isempty(ProbFormTemp)
            [AMatTemp,BrhsTemp] = OptimiFull(ProbFormTemp);
            
            if ~isempty(AMatTemp)
                if (size(AMatTemp,2) ~= nD)
                    warning(['nearoptplotmo2: # columns of ProbForm data inconsistent with other data passed. Reverting to default of empty constraints'])
                else
                    AMat = AMatTemp;
                    Brhs = BrhsTemp;                    
                    cFunc=cFuncTemp;
                    blValidLinProgramData = 1;
                    
                    %Set empty fields of ProbForm to default values of []
                    %for convienence later in retrieving and working with
                    %values
                    if ~isfield(ProbFormTemp,'Aineq')
                        ProbFormTemp.Aineq = [];
                        ProbFormTemp.bineq = [];
                    end
                    if ~isfield(ProbFormTemp,'Aeq')
                        ProbFormTemp.Aeq = [];
                        ProbFormTemp.beq = [];
                    end
                    if ~isfield(ProbFormTemp,'lb')
                        ProbFormTemp.lb = [];
                    end
                    if ~isfield(ProbFormTemp,'ub')
                        ProbFormTemp.ub = [];
                    end
                    
                    %Set default fields for solver and options if not
                    %defined
                    if ~isfield(ProbFormTemp,'solver')
                        ProbFormTemp.solver = 'linprog';
                    end
                    if ~isfield(ProbFormTemp,'options')
                        ProbFormTemp.options = struct('maxiter',15000,'Display', 'off','Algorithm','interior-point');
                    end
                   
                    ProbForm = ProbFormTemp;
                end
            end
        end
        
         %Further error testing on parameters needed to generate alternatives       
         if (GenMethod==2) && (isempty(AMat) || isempty(Brhs) || isempty(NearOptConstraint) || isempty(OptSolRow))
            warning('nearoptplotmo2: To generate alternatives by Matlab linear programming, must define AMat, Brhs, NearOptConstraint, and OptSolRow inputs. Setting GenerateMethod engine to Data and other inputs to empty.');
            GenerateMethod = 1;
            AMat = []; Brhs=[]; NearOptConstraint=[]; OptSolRow=[];
         end         
         
         if (GenMethod==3) && (isempty(sGamsFile) || ~exist(sGamsFile,'file'))
            warning(['nearoptplotmo2: ',sGamsFile, ' GAMS file does not exist. Setting to empy. To generate alternatives by GAMS, enter a valid file name on the Interact Tab.']);
            sGamsFile = '';            
         end           
         
         %If showing group labels table, need to have the table of
         %mActCat to show
         if lShowGroupLabels && strcmpi(lShowGroupLabelsEnable,'off')
             warning('nearoptplotmo2: Need to specify the mActCat of decision axes groupings to show the grouping table. This feature will be disabled.')
             lShowGroupLabels = 0;
         end
    end
    
    %Error checking on vGroup
    
    %if vGroup is inumeric, convert it to a cell array of strings
    if isnumeric(vGroup)
        vGroup = arrayfun(@num2str,vGroup, 'unif', 0);
    end
    
    %Error checking and default settings on Current Record and iOpt      
       if ~strcmpi(iOpt,'0')
          %Set current record to the record of the highlight group
          lCurrRecord = find(ismember(vGroup,iOpt)~=0,1,'first');
        else
          if lCurrRecord < 1
              lCurrRecord = 1;
          elseif lCurrRecord > m
              lCurrRecord = m;
          end
       end
    
       lHighRecord = lCurrRecord; %Also save as the highlighted record
    
        
    %Determine the number of groups.
    [vUnique vIC vID] = unique(vGroup);
    nU = max(size(vUnique));
    vShowGroup = ones(nU,1);
    vGroupThick = ones(nU,1);
    
    if ~isempty(mGroupData)
       %Compare vUnique to vGroupOrder 
       IsMemberInd = ismember(mGroupData(:,1),vUnique);
       if (sum(IsMemberInd)==length(vUnique)) && (length(vUnique)==size(mGroupData,1)) && (size(mGroupData,2)==3) %Check correct # rows and # cols
           vUnique = mGroupData(:,1);
           vShowGroup = (cell2mat(mGroupData(:,2))>0);  
           vGroupThick = (cell2mat(mGroupData(:,3))); 
       else
           warning('mGroupData has different members than vGroup or wrong number of columns. Reverting to default alpha-numeric sorting of groups and showing all groups')
       end
    end
    
    lOptGroupCurr = find(ismember(vUnique,iOpt)~=0,1,'first');
    if isempty(lOptGroupCurr)
        lOptGroupCurr = 0;
    end
    
    %Error Check the color parameters
    colorerrormsg = '';
    
    if ~isempty(mColors) && ((size(mColors,2) < 1 +(nO>0)) || (size(mColors,3)<3) || (size(mColors,4) ~= 3))
        colorerrormsg = sprintf('mColors is of incorrect size. Must be %d (# Groups) by %d by 3 (color depth) by 3 (R-G-B) matrix. Using default colors.',nU,1+(nO>0));
        mColors=[];
    end

    if ~isempty(mColors) 
        vUniqueGroupColors = unique(mColors(:,1));

        if nU > length(vUniqueGroupColors)
            colorerrormsg = sprintf('Number of groups exceeds number of group colors provided. Will reuse colors for group #s %d and above.',nU+1);
        end
    end

    if ~strcmpi(colorerrormsg,'')
        warning(['nearoptplotmo2: mColors. ',colorerrormsg]);
    end

    if isempty(mColors)
        %Build the default color matrix           
        mColors = zeros(7,2,3,3);
        %Group 1 gets Green To Magenta
        GtoM = OSUColorRamps('GreenToMagenta16Step');
        %Decisions -> the greens
        mColors(1,1,1:3,:) = flipud(GtoM([1 3 5],:));
        %Objectives -> the magentas
        mColors(1,2,1:3,:) = GtoM([11 13 15],:);

        %Pull groups 2-5 from the stepped sequential scheme
        StepSeq = flipud(OSUColorRamps('SteppedSequential'));           
        ssRowInds = [2 4 5 17 19 20 7 9 10 22 24 25]';

        for i=2:5
            mColors(i,1,[1:3],:) = StepSeq(ssRowInds(3*(i-2)+1:3*(i-1)),:);
        end
        %Pull groups 6 and 7 from  categorial 12 step       
        CatStep = ExpandColorRamp(1,OSUColorRamps('Categorical12Step'),1,0); %expand the ramp to get 3 color hues per group             
        for i=6:7
            mColors(i,1,[1:3]',:) = StepSeq(3*(i-6)+[1:3],:);      
        end
        %copy over dec to objs for groups #2 up
        mColors([2:7],2,:,:) = mColors([2:7],1,:,:);
    end
    
       %If there are objective function axes, we can enable the option
       %to color them in a contrasting color
       if nO>0
           sEnableObjColors = 'on';
           lObjColorsVis = lShowObjsDiffColor;
       else
           sEnableObjColors = 'off';
           lObjColorsVis = 0; %set as unchecked but will pass the original setting (lShowObjsDifferentColor) back out
       end
       
    % Error checking on the colors for axis labels and the grouping table
    [vUniques, iC, uIndex, iCBeg, iVBeg] = GroupOrderAppear(mActCat(:,1));
    nCats = length(vUniques);
    nF = size(mActCat,2);
    
    if ~isempty(mColorsAxisLabels)
        nSs = size(mColorsAxisLabels);
        if (length(nSs)~=3) || (nSs(1) ~= nCats) || (nSs(2) ~= nF) || (nSs(3) ~= 3) 
            colorerrormsg = sprintf('mColorsAxisLabels is of incorrect size. Must be %d x %d x 3 matrix where %d=number of categories in first field of mActAct, %d=number of fields in mActAct, and 3=RGB values. Using default colors.',nCats,nF,nCats,nF);
            warning(colorerrormsg)
            mColorsAxisLabels=[];
        end
    end

    if isempty(mColorsAxisLabels)
        %Draw from OSU Stepped Sequential Scheme
        mStepSeq = OSUColorRamps('SteppedSequential');
        %Transform 2D in into 3D
        mColorFull = zeros(nCats,nF,3);
        cCats = repmat([2 5 1 3 4],1,ceil(nCats/5)); %orange-purple-red-green-blue family of stepped sequential, wrap to repeat
        for i=1:nCats
           if nF <=5
                mColorFull(i,:,:) = mStepSeq(5*(cCats(i)-1)+[1:nF],:); 
           else
               mColorFull(i,:,:) = ExpandColorRamp(nF-2,mStepSeq(5*(cCats(i)-1)+[1:nF],:),4,0);
           end
        end
    else
        mColorFull = mColorsAxisLabels;  
    end
    
    mColorChoices = mColorFull(:,nF,:);
    
    %Error checking on colors for left and right ticks and axis labels
    if ~isempty(mColorYScales)
        nCTs = size(mColorYScales);
        if (length(nCTs)~=2) || (nCTs(1) ~= 2) || (nCTs(2) ~= 3) 
            colorerrormsg = sprintf('mColorYScales is of incorrect size. Must be 2 x 3 matrix where row 1 is for left ticks, row 2 is for right ticks, 3=RGB values. Using default colors.');
            warning(colorerrormsg)
            mColorYScales=[];
        end
    end   
    
    if isempty(mColorYScales)
       mColorYScales = reshape(mColors(1,1:2,3,:),2,3); 
    end
    
    % Error checking on the precision -- round to nearest integer
    Precision = round(Precision);
    
    % Round the mResults to the specified precision
    mResults = round(mResults*10^Precision)/10^Precision;
    
    % update vStep if no value passed
    if isnan(vStep)
        vStep = (max(mResults,[],1)-min(mResults,[],1))/20;
    end
    
    %% Start building the plot
     
    %Get the Figure Handle
    if (hWind>0) && ishandle(hWind)
        Figure2 = hWind;
        clf(Figure2);
        %sprintf('Cleared figured')
    else
        Figure2=figure('units','normalized','OuterPosition',[0.02 0.04 .96 .90]);
        %sprintf('Maximize figure')
    end
    
    hWindReturn = Figure2;
       
    lLineWidth = 1.5;
    
    %Determine the size of the figure components (Parallel Plot and Control
    %Panel)
    %Shift to pixel units to get current figure size
    cUnits = get(hWindReturn,'Units');
    set(hWindReturn,'Units','characters');  %set(hWindReturn,'Units','pixel');
    vSizePx = get(hWindReturn,'OuterPosition');
    set(hWindReturn,'Units',cUnits);
    
    lMarginPanel = 0.1; % margin between parallel coordinate plot and control panel in normalized units
    rMarginPanel = 0.01; % margin between control panel and figure right edge in normalized units
    xPanelWidthNorm = PanelWidth/vSizePx(3); %width of control panel in normalized units
    
    %Covert the plot position to local variables
    xLeft = PlotPosition(1);
    xWidth = 1 - (xLeft + lMarginPanel + xPanelWidthNorm + rMarginPanel); %PlotPosition(3);
    yBottom = PlotPosition(2);
    yHeight = PlotPosition(4);
    PlotPosition(3) = xWidth;
    
    plot2 = axes('Parent',Figure2,'FontSize',FontSize-4,'Units','normalized','Position',[xLeft yBottom xWidth yHeight],'YColor',mColorYScales(2,:),'box','on');
    
    hold on
    
    %get(plot2,'Position')
    
    %fprintf('Initial Maxes / Mins\n');
    vMaxResult = max(mResults,[],1);
    vMinResult = min(mResults,[],1);
    
    %% Rescale/transform the vertical axis if needed
    % We need to find a set of plot units that allow us to show the data in
    % each column (potentially of different magnitudes and ranges) on a common plot
    
    %Start off by specifying the lower and upper y plot limits
    UseStandarize = 'off';
    vMaxPlot = vMaxResult;
    vMinPlot = vMinResult;
    
    switch (AxisScales)
        case 'none'
            BaseAxis = [1 1];
            ymaxdec = max(vMaxResult); 
            ymin = min(vMinResult);
            
        case 'standardize'
            UseStandarize = 'on'; 
            %calc the max deviation from the mean for each axis
            vMeans = mean(mResults);
            vDev = std(mResults);
            vMaxDev = zeros(1,n);
            for i=1:n
                vMaxDev(i) = max(abs(mResults(:,i)-vMeans(i))/vDev(i));
            end
            
            ymaxdec = max(vMaxDev);
            ymin = -ymaxdec;
        
        case 'normalize'
            ymaxdec = 1;
            ymin = 0;
            
        case 'auto'
            [ymaxdec, iDex] = max(vMaxResult(nO+1:n)); 
            iObjDex = nO;
            if nO>0
                [yobjmax, iObjDex] = max(vMaxResult(1:nO));
            end
            BaseAxis = [iObjDex iDex];
            ymin=0;
            
        case 'custom'
            ymaxdec = mLims(2,nO+BaseAxis(2));
            ymin = mLims(1,nO+BaseAxis(2));
            
    end
    
    vLims = [1-YAxisMargins(1) n+YAxisMargins(2)];
    %set the limits and get the ticks   
    if 1
        set(plot2,'xLim',vLims);
        %dTemp = plot([ymin ymaxdec]);
        %yLimsStart = get(plot2,'yLim');
        %dticks = [yLimsStart(1):(yLimsStart(2)-yLimsStart(1))/(NumTicks-1):yLimsStart(2)];
        %delete(dTemp);
    else %old approach
        dticks = [ymin:(ymaxdec-ymin)/(NumTicks(1)-1):ymaxdec];
        set(plot2,'yLim',[ymin ymaxdec],'xLim',vLims,'ytick',dticks,'yticklabel',dticks);
    end
    %dticks = str2num(get(plot2,'yticklabel'));   
    
    %scale the objective and decision variable axes
    %vObjScale is a multiplier by which to multiply the axis to get it to
    %fit within the plot area
    vAxisScaleAll = ones(1,nO+nD);
    
    if nO>0
        vAxisLabelsAll = {vObjLabels{:} vXLabels{:}};
    else
        vAxisLabelsAll = vXLabels;
    end
    
    if nO==0
        yAxisLabels{1} = yAxisLabels{2};
    end
    
    AxisLabelsNew = yAxisLabels;
    
    mPlot = mResults;
    ymaxobj = 0;
    
    adjust = 1;
       
    %now make the transformations
    if strcmpi(AxisScales,'normalize')
        for i=1:n
            mPlot(:,i) = (mResults(:,i) - vMinResult(i))/(vMaxResult(i)- vMinResult(i));
        end
        mResults;
        yLimsSecond = [0 1];
        
     elseif strcmpi(AxisScales,'standardize')
        for i=1:n
            mPlot(:,i) = (mResults(:,i) - vMeans(i))/(vDev(i));
        end
        
        mResults; 
        UseStandarize = 'off';
        
     elseif strcmpi(AxisScales,'auto') || strcmpi(AxisScales,'custom')           
           BaseAxis;
           
           vAdjust = zeros(n,1);
           TransOffsets = zeros(n,1);
           minMults = zeros(1,2);
           vMults = zeros(1,n);
           
           if strcmpi(AxisScales,'auto')
                vMaxCompare = vMaxResult;
                mPlotTemp = mPlot;

                %Now build the entire set of transformed objectives and decision axes
                [mUseResAll, vMultAll, BaseAll] = ScaleMultipleAxes(mPlotTemp);
                
                mPlot = mUseResAll;
                BaseAll;
                
                MaxAllVal = max(mUseResAll(:,BaseAll),[],1);
                %Concatenate all the multipliers
                vMults = vMultAll;
                %vMults = [vMultObj vMultAll(nO+1:n)];
                [minMults(2), minDC] = min(vMults(nO+1:n));
                
                %MaxToUse = MaxAllVal;
                %MinToUse = min(min(mUseResAll));
                %Construct the ticks

            else
               vMaxCompare = mLims(2,:);
               mLimsTransformed = mLims;
               mPlotTemp = mLims;
               
               %First pass 
               for i=1:n            
                   CompAxis = nO+BaseAxis(2);
                   
                   TransSlope(i) = (mLims(2,CompAxis)-mLims(1,CompAxis))/(mLims(2,i)-mLims(1,i));
                   TransOffsets(i) = +mLims(1,CompAxis) - TransSlope(i)*mLims(1,i);
                   mPlot(:,i) = TransSlope(i)*mResults(:,i)+TransOffsets(i);
                   
                   if i<=nO
                       CompAxis = BaseAxis(1);
                   else
                       CompAxis = BaseAxis(2);
                   end
                   
                   vAdjust(i) = (mLims(2,i)/mLims(2,CompAxis));
                   vMultAll(i) = floor(log(mLims(2,i))/log(10));
               end
               
               vAdjust';
               TransSlope;
               TransOffsets;
               %lticks;
               %vMultAll = log(vAdjust')/log(10)
               %objticks = objticks*10^min(vMultAll);
               %vMultAll = vMultAll-min(vMultAll);
              
               vMults = vMultAll-1;
               
               for k=1:2
                    minMults(k)=vMults(nO*(k-1)+BaseAxis(k));
               end
               
               if nO==0
                   minMults(1) = minMults(2);
               end
                 
           end
    end
    
    %Calculate the transformation from plot units back to original data units
    %matrix to store directions for how to transform plot units back to the data units. 
    %Row 1 is the slope, Row 2 the intercept, i.e., DataUnits = Slope*(PlotUnits) + Intercept
    mTransformToOrig = zeros(2,n);
    
    if strcmpi(AxisScales,'custom')
        %We already calculated the coeeficients
        mTransformToOrig(1,:) = 1./TransSlope(:);                   
        mTransformToOrig(2,:) = -TransOffsets(:)./TransSlope(:);
        
        ymaxPU = ConvertPlot2DataUnits(mLims(2,nO+BaseAxis(2)),mTransformToOrig(:,nO+BaseAxis(2)),1);
        yminPU = ConvertPlot2DataUnits(mLims(1,nO+BaseAxis(2)),mTransformToOrig(:,nO+BaseAxis(2)),1);
        yLimsSecond = [yminPU ymaxPU];
    else    
        %Use linear regression against the mResults and mPlot data in each
        %column
        for i=1:n
            mTransformToOrig(:,i) = polyfit(mPlot(:,i),mResults(:,i),1)';
        end
    end
    mTransformToOrig;    
    
    %Revisit the scaling now in plot units   
    %Temportary plot in plot units to get the auto Y limits
    dTemp = parallelcoords(mPlot);
    if ~strcmpi(AxisScales,'custom') && ~strcmpi(AxisScales,'normalize')
        %allow Matlab to determine the y axis limits            
        yLimsSecond = get(plot2,'yLim');
    end
    yLimsSecond;
    ymax = yLimsSecond(2);
    ymin = yLimsSecond(1);
    %Create a set of ticks along these limits
    dticks = [yLimsSecond(1):(yLimsSecond(2)-yLimsSecond(1))/(NumTicks(1)-1):yLimsSecond(2)];
    dticksr = [yLimsSecond(1):(yLimsSecond(2)-yLimsSecond(1))/(NumTicks(2)-1):yLimsSecond(2)];
    set(plot2,'ytick',dticks);
    delete(dTemp);
    
    if strcmpi(AxisScales,'auto') || strcmpi(AxisScales,'custom')
        %May need to revisit scale labels/ticks on the left and right scales
        vMults;
        minMults;
        lticks=dticks; %left for objective scale
        rticks=dticksr; %right for decision scale
           
        if strcmpi(AxisScales,'auto')

            rticks = rticks*10^(minMults(2));
            if nO>0
                [minMults(1), minOC] = min(vMults(1:nO));
                %Construct the ticks from the min column
                lticks = lticks*10^(minMults(1));                   
            else
                minMults(1) = minMults(2);
                minOC = minDC;
                lticks = rticks;
            end
        else
           for i=1:2
               %currticks = [mLims(1,i):(mLims(2,i)-mLims(1,i))/(length(dticks)-1):mLims(2,i)];
               currticks = [mLims(1,i):(mLims(2,i)-mLims(1,i))/(NumTicks(i)-1):mLims(2,i)];
               
               if i==BaseAxis(1)
                   lticks = currticks;
               elseif i==nO+BaseAxis(2)
                   rticks=currticks;
               end
           end
        end
        
           %mResults;
           %mPlot;

           %dticks;
           %lticks;
           %rticks;
                     
            %work on the scale and label for the left (objective) and right (decision) axis labels in turn                                              
            for k=1:2
                if (k==1)
                    currticks = lticks;
                    currplot = plot2;
                    cStart = 1; cEnd = nO; %start, end columns
                else
                    currticks = rticks;
                    %currplot = ax2;
                    cStart = nO+1; cEnd = n; %start, end columns
                end
                
                if nO==0
                    cStart = 1; cEnd = n;
                end
                
                if minMults(k)>=TickMag; %5
                    %move the surplus zeros into the axis label
                    %first read the label
                    if isempty(findstr(yAxisLabels{k},'('))
                        s1=yAxisLabels{k};
                        s2=')';
                    else
                        [s1 s2] = strread(yAxisLabels{k},'%s %s','delimiter','(');
                        s2 = [' ',s2];
                    end
                    %reformat the label with adjustment factor added

                    baseP = minMults(k);
                    AxisLabelsNew{k} = sprintf('%s (10^{%d}%s',char(s1),baseP,char(s2));
                    currticks=currticks/10^minMults(k);
                %elseif minMult>=2
                    %move the surplus zeros into the ticks
                    %ObjTicksFinal = ThousandSep(round(objticks));
                %     objticks = objticks*10^minMult;
                end

                if (minMults(k)~=0)
                    %Update the ticks
                     TicksFinal = ThousandSep(currticks);    
                     if k==1
                        set(currplot,'yticklabel',TicksFinal); 
                     else
                         %right axis doesn't yet exist, save ticks
                         rticks = TicksFinal;
                     end

                    %Transform the plot data
                    if (nO>0) || ((nO==0) && (k==2))
                        vMults(cStart:cEnd) = vMults(cStart:cEnd)-minMults(k);
                    end
                    %vMultObj = vMultObj-minMult;
                    %mPlot(:,cStart:cEnd) = mPlot(:,cStart:cEnd)/(10^minMults(k));
                end                      
            end

            %vMults = [vMultObj vMultAll(nO+1:n)-minMultDec]; %ignore the multipliers on the objective function axes since we already dealth with them.
            vMults;
            minMults;
            %adjust the obj and decision axes labels if needed based on the
            %scale multiplier
            for i=1:n
                if (vMults(i)~=0) 
                      %first read the label
                      if isempty(findstr(vAxisLabelsAll{i},'('))
                         s1=vAxisLabelsAll{i};
                         s2='';
                     else 
                         [s1, s2] = strread(vAxisLabelsAll{i},'%s %s','delimiter','(');
                         s2 = strcat(' (',s2);
                     end
                    %reformat the label with adjustment factor added
                    if TransOffsets(i) == 0
                        if (vMults(i)>3) || (vMults(i) < 0)
                            vAxisLabelsAll{i} = sprintf('%s [10^{%d}%s]',char(s1),vMults(i),char(s2));
                        else
                           vAxisLabelsAll{i} = sprintf('%s [x%.0f%s]',char(s1),10^vMults(i),char(s2)); 
                        end
                    else
                       vAxisLabelsAll{i} = sprintf('%s [x%.0f + %.0f]',char(s1),10^vMults(i),TransOffsets(i),char(s2));
                    end

                end
            end      
           
    end

   % if strcmpi(AxisScales,'normalize')
        %vMaxPlot(:) = 1;
        %vMinPlot(:) = 0;  
    %else
        vMaxPlot = max(mPlot,[],1);
        vMinPlot = min(mPlot,[],1);
   % end
    
    mPlot(:,1:nO);
    
    %Create a vector to store vFixedVals in Plot Units
    vFixedValsPlotUnits = vFixedVals;
    %correct vFixedVals if outside of range. Also 
    for i=1:n
        if vFixedVals(i) > vMaxResult(i)
            vFixedVals(i) = vMaxResult(i);
        end
        if vFixedVals(i) < vMinResult(i)
            vFixedVals(i) = vMinResult(i);
        end
        
        vFixedValsPlotUnits(i) = ConvertPlot2DataUnits(vFixedVals(i),mTransformToOrig(:,i),1);
        
        if vFixedValsPlotUnits(i) > vMaxPlot(i)
            vFixedValsPlotUnits(i) = vMaxPlot(i);
        end
        if vFixedValsPlotUnits(i) < vMinPlot(i)
            vFixedValsPlotUnits(i) = vMinPlot(i);
        end
    end
       
    hTexts = zeros(1,n);
    sLabel = vObjLabels;
    
    %% Start Working on the Rows/lines
    % Rows are members of a group (as specified by vGroup
    % and each group is plotted in a different color
        
    lMainClass = 1;
    mTextPos = zeros(n,4);
    chHeight = 1.3; %in characters
    
    %display the axis lines and labels
    for i=1:n    
           %['i: ', num2str(i),'; nO: ', num2str(nO)]
           
            h3 = plot([i i],[ymin ymax]);
            if i>nO %decisions
                lLineColor = [0.8 0.8 0.8];
                lMainClass = iVBeg(i-nO);
                vColor = mColorChoices(iVBeg(i-nO),:);
            else
                lLineColor = [0.6 0.6 0.6];
                vColor = [.737 0 0.737];
            end
            set(h3,'linewidth',0.5,'color',lLineColor);
            hTexts(i) = text(i,0,vAxisLabelsAll{i},'fontsize',FontSize-6,'HorizontalAlignment', 'right','color',vColor,'VerticalAlignment','middle','Rotation',90);
            %set(hTexts(i),'EdgeColor',[0 0 0]); %for help in placing group labels box below
            %ymin-6*(ymax-ymin)/100
           %Shift the text down belwo the checkbox
            strCurrUnits = get(hTexts(i),'Units');
            set(hTexts(i),'Units','characters');
            hTextPos = get(hTexts(i),'Position');
            hTextPos(2) = hTextPos(2) - chHeight-1;
            set(hTexts(i),'Position',hTextPos);
            set(hTexts(i),'Units','normalized');
            mTextPos(i,:) = get(hTexts(i),'Extent');
            set(hTexts(i),'Units',strCurrUnits);
    end
    
     
     %mTextPos
    %Alternative way of labeling decision variable axes based on mActCat
    %old data units
    %[hTextGroup, rDex, hGLines] = GenerateGroupLabels(mActCat,1,0,[1  ymin-6*(ymax-ymin)/100 3.3*(ymax-ymin)/100 0 4*(ymax-ymin)/100],'data',FontSize-2);
    
    yBottomTable = 0.025;
    mMaxExtent = min(mTextPos(:,2));
    yTableHeight = max([yBottom-abs(yHeight*mMaxExtent)-yBottomTable-0.005; 0.02*(sum(FontSize-4-12+[1:size(mActCat,2)-1])); 0.001]);
    %position the table below the largest axis text label
    vTablePosition = [xLeft+(nO-1+0.5-.0025)/(n-0.5)*xWidth yBottomTable xWidth*(n-nO+0.25)/(n-0.5) yTableHeight];
    %[yBottomTable yTableHeight]
    %vTablePosition = [xLeft+(YAxisMargins(1)+nO-1+0.5-.0025)/(n-1+sum(YAxisMargins))*xWidth yBottomTable xWidth*(nD)/(n+sum(YAxisMargins)) yTableHeight];
    lFieldExclude = 0;
    
    if lFieldExclude == 0
        lNumFields = nF;
    else
        lNumFields = lFieldExclude-1;
    end
    
    [hTabComponent hTabContainer]=GenerateGroupLabelsTable(hWindReturn,mActCat,lFieldExclude,vTablePosition,'normalized',mColorFull,FontSize-3);
    
    hTextGroup = []; rDex = []; hGLines = [];
    
    
    [nActR nActC] = size(mActCat);
    yFirstHeight = 0.25*yBottom;
  % [hTextGroup, rDex, hGLines] = GenerateGroupLabels(mActCat,1,0,[xLeft+nO/(n+0.5)*xWidth yBottom-.05 (yBottom-yFirstHeight)/(nActC-1) yFirstHeight 0.025 xWidth],'normalized',FontSize-2);
    
    hTGSize = max(size(hTextGroup));
    for i=1:hTGSize
       cPos = get(hTextGroup(i),'Position');
       %sprintf('%d %s %d %.2f\n',i,get(hTextGroup(i),'string'),cPos(1),cPos(2)) 
    end
       
    %% Plot the groups

    if 1 %nU > 1
        %overplot additional groups in a separate color/formatting
        %first pull out the group
               
        mColorsDecsGroup = mColors(:,1,1,:);
        mColorsObjsGroup = mColors(:,1+(nO>0),1,:);
        
        vLineStyle = {'-' '-' '-' '-' '-' '-' '-'};
        %vLineStyle = {'-' '-' '-.'};

        hPCGroup = cell(nU,2); %Create a cell array of handles to the plots to be created nU groups by (decision variables and objectives) 
       
        for i=1:nU %Loop through the groups            
            if iscell(vUnique(i))
                mGroup = mPlot(strcmpi(vGroup,vUnique(i)),:);
            else
                mGroup = mPlot(vGroup==vUnique(i),:);
            end
            
            [lIndToUse,lLineWidth,mColorsDecs,mColorsObjs,strVis] = GetGroupRenderings(i,[]);
          
            hPCGroup{i,1}= parallelcoords(mGroup,'color',mColorsDecs,'LineStyle',vLineStyle{lIndToUse},'linewidth',lLineWidth,'Standardize',UseStandarize,'Visible',strVis);
            if (nO>0) %overplot the objective axes in a different color
                strVisObj = Boolean2Enabled(vShowGroup(i) && lObjColorsVis);
                hPCGroup{i,2}= parallelcoords(mGroup(:,1:(nO+1)),'color',mColorsObjs,'LineStyle',vLineStyle{lIndToUse},'linewidth',lLineWidth,'Standardize',UseStandarize,'Visible',strVisObj);
            end          
        end
    else
        %no grouping -- plot all rows in the same color
        
        lLineWidth=vGroupThick(1);
        %vGroup = ones(1,m);
        h2 = parallelcoords(mPlot,'color',[0.526 1 0.526],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize);
        if (nO>0) %overplot the objective axes in a different color
            strVisObj = Boolean2Enabled(lObjColorsVis);
            h3 = parallelcoords(mPlot(:,1:(nO+1)),'color',[1 0.526 1],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize,'Visible',strVisObj);
        end
    end
    
    %% Cleanup the plot
    
    set(plot2,'yLim',[ymin ymax],'Xticklabel',[],'Xtick',[],'fontsize',FontSize-4);

    %print the title + axis labels
    mytitle = sprintf('Near optimal region');
    %title(mytitle,'FontSize',FontSize);
    ylabel(plot2, AxisLabelsNew{1},'fontsize',FontSize);
    
    %duplicate axes to plot decision variable values at right
    ax2 = axes('Position',get(plot2,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'YColor',mColorYScales(1,:));
   
    set(ax2,'yLim',[ymin ymax],'xLim',vLims,'Xticklabel',[],'Xtick',[],'fontsize',FontSize-4);
    RightTicks = str2num(get(ax2,'YTickLabel'));
    class(RightTicks);
    RightTicks;
    ylabel(ax2,AxisLabelsNew{2},'fontsize',FontSize);

    
    %update the tick labels on the left and right scales
    
    
    if strcmpi(AxisScales,'normalize')   
        nLeftTicks = num2cell(dticks);
        nRightTicks = nLeftTicks;
        nLeftTicks{1} = 'Min. - 0';
        nLeftTicks{end} = 'Max. - 1.0';
        %nRightTicks = num2cell(RightTicks);
        nRightTicks{1} = '0 - Min.';
        nRightTicks{end} = '1.0 - Max.';
    elseif strcmpi(AxisScales,'standardize')
        nLeftTicks = num2cell(RightTicks);

        nTicks = length(nLeftTicks);
        nRightTicks = nLeftTicks;
        for i=1:nTicks
           if floor(RightTicks(i))==RightTicks(i)
               strWord = 'S. Dev.';
               strAdd = ' ';
               if RightTicks(i)==0
                   strWord = 'Mean';
                   strAdd = ' - ';
               end
               nRightTicks{i} = [num2str(RightTicks(i)), strAdd, strWord];
               nLeftTicks{i} = [strWord,strAdd, num2str(RightTicks(i))];
           end
        end
    elseif strcmpi(AxisScales,'none')
        %nRightTicks = ThousandSep(get(ax2,'ytick'));
        nLeftTicks = ThousandSep(get(plot2,'ytick'));
        nRightTicks = nLeftTicks;
    elseif strcmpi(AxisScales,'custom') || strcmpi(AxisScales,'auto')
        nLeftTicks=get(plot2,'yTickLabel');
        nRightTicks = rticks';
    else
        nRightTicks = ThousandSep(RightTicks);
        nLeftTicks=get(plot2,'yTickLabel');
    end
    %plot last one as objective function
    nLeftTicks;
    nRightTicks;
   
    %Calc and set new ticks
    lYLim = get(plot2,'yLim');
    rYLim = get(ax2,'yLim');
    
    [length(nLeftTicks) length(nRightTicks)];
    lTickVals = [lYLim(1):(lYLim(2)-lYLim(1))/(length(nLeftTicks)-1):lYLim(2)];
    rTickVals = [rYLim(1):(rYLim(2)-rYLim(1))/(length(nRightTicks)-1):rYLim(2)];
    
    
    set(plot2,'yTick',lTickVals,'yTickLabel',nLeftTicks);
    set(ax2,'yTick',rTickVals,'yTickLabel',nRightTicks);
      
    % Determine the ticks on the right axis scale
    TicksRight = get(ax2,'ytick');
    TickSize = TicksRight(2)-TicksRight(1);
    
    %Rebuild the variable argument list to pass along through the controls
    
%    if strcmpi(AxisScales,'custom')
        vararg_curr = {'hWind' hWindReturn 'Tolerance' Tolerance 'AllowableDeviation' AllowableDeviation ...
                'ErrorResid' ErrorResid 'NumSamples' NumSamples 'ObjSamplesPercent' ObjSamplesPercent  ...
                'vFixed' vFixed 'Precision' Precision 'vFixedVals' vFixedVals 'vStep' vStep 'FontSize' FontSize 'NumTicks' NumTicks 'TickMag' TickMag 'mActCat' mActCat ...
                'vGroup' vGroup 'mGroupData' mGroupData 'GroupToHighlight' iOpt 'mColors' mColors 'mHighlightColor' mHighlightColor ...
                'ProbForm' ProbForm 'cFunc' cFunc 'PlotPosition' PlotPosition 'PanelWidth' PanelWidth 'vObjLabels' vObjLabels ...
                'vXLabels' vXLabels 'vXLabelsShort' vXLabelsShort 'yAxisLabels' yAxisLabels 'AxisScales' AxisScales 'BaseAxis' BaseAxis 'mLims' mLims ...
                'TickSize' TickSize 'sGamsFile' sGamsFile 'StartTab' StartTab, 'GenerateType' GenType 'GenerateMethod' GenMethod ...
                'HideSliders' HideSliders 'HideCheckboxes' HideCheckboxes 'NearOptConstraint' NearOptConstraint 'OptSolRow' OptSolRow ...
                'mColorsAxisLabels' mColorsAxisLabels 'mColorYScales' mColorYScales 'YAxisMargins' YAxisMargins};

         %Also set app variables for the functon inputs and computations
         setappdata(hWindReturn,'varargs',vararg_curr);
         setappdata(hWindReturn,'mConvert',mTransformToOrig);
         setappdata(hWindReturn,'mDecs',mDecisions);
         setappdata(hWindReturn,'mObjs',mObjs);

    %Check All Axes Checkbox
    if sum(vFixed)==n
        AllValue = 1;
    else
        AllValue = 0;
    end    

%% Add controls to the plot
% There are 5 basic elements
% A) Menu items to the right of existing menus to allow (i) saving data and
%       loading results as well as (ii) toggling/hiding the other controls,
%       and (iii) updating the model formulation
% B,C,D) 3 elements are frames for the Interact, Color Ramp, and Display tabs that have lots
%     of controls
% E) The final elements comprise the siders and axis check boxes on the plot
%     itself         

   %Data menu with save data and load results options
   dMenu = uimenu('Label','Plot Data');
      uimenu(dMenu,'Label','Save data','Callback',{@SaveFig,hWindReturn},'Accelerator','D');
      uimLoadData = uimenu(dMenu,'Label','Load data...','Enable','on','Callback',{@LoadData,hWindReturn},'Accelerator','L');
      uimPrintToPDF = uimenu(dMenu,'Label','Print to PDF...','Callback',{@PrintToPDF,hWindReturn,mResults,nO});
   
   %Control menu with options to show/hide/set all the controls, sliders,
   %checkboxes, plot layers showing various categories of alternatives, coloring for values on objective function axes, and
   %table below the axes labels to group the axes according to mActCat.
   %Additionally, give the user a second place to go directly to the
   %Interact, Color Ramp, and Display tabs
   cMenu = uimenu('Label','Controls');
       uimShowAll = uimenu(cMenu,'Label','Hide all controls','Checked','off');
       uimSliders = uimenu(cMenu,'Label','Sliders','separator','on');
          uimHideSliders = uimenu(uimSliders,'Label','Hide sliders','Checked',Boolean2Enabled(HideSliders));
          uimSetBeyond = uimenu(uimSliders,'Label','Set beyond extents','Enable','off');
       uimCheckboxes = uimenu(cMenu,'Label','Axes checkboxes');
          uimHideChecks = uimenu(uimCheckboxes,'Label','Hide checkboxes'); %'Checked',Boolean2Enabled(HideCheckboxes)); Don't set here. Later will trigger with callback
          uimCheckAll = uimenu(uimCheckboxes,'Label','Check all','Checked',Boolean2Enabled(AllValue));
       uimShowAlts = uimenu(cMenu,'Label','Alternatives'); 
          uimShowCurrRec = uimenu(uimShowAlts,'Label','Show current','Checked',Boolean2Enabled(lShowCurrRecord)); 
          uimShowFiltered = uimenu(uimShowAlts,'Label','Show filtered','Checked',Boolean2Enabled(lShowFiltered));       
       uimShowObjsDiffColor = uimenu(cMenu,'Label','Contrast colors for objectives','enable',sEnableObjColors,'checked',Boolean2Enabled(lObjColorsVis),'Separator','on');
       uimGroupLabels = uimenu(cMenu,'Label','Show group labels','enable',lShowGroupLabelsEnable); %Don't set check value here; only later should lShowGroupLabels=1 and we trigger the event callback
       uimShowInset =   uimenu(cMenu,'Label','Show inset plot','enable',Boolean2Enabled(nO==2)); %Don't set check value here; only later should ShowInset==1 and we trigger the event callback
       uimShowLegend =  uimenu(cMenu,'Label','Show legend'); 
       uimLegendOptions = zeros(numel(cLegendOptions),1);
       for i=1:numel(cLegendOptions) %Mutually exclusive
          uimLegendOptions(i) = uimenu(uimShowLegend,'Label',cLegendOptions{i},'Checked',Boolean2Enabled(strcmpi(cLegendOptions{i},ShowLegend)),'Callback',{@UpdateLegendMenu,i,numel(cLegendOptions)}); 
       end

       %Menu items to go directly to the tabs
       uimTab = zeros(length(ButtonText),1);
       for i=1:length(ButtonText)
           if i==1
               uimTab(i) = uimenu(cMenu,'Label',ButtonText{i},'Separator','on');          
           else
               uimTab(i) = uimenu(cMenu,'Label',ButtonText{i});
           end
       end

       uimRefresh = uimenu(cMenu,'Label','Refresh display','Separator','on','Callback',{@ReorderGroups,mResults,nO,hWindReturn});
       
   %Update model formulation tab to allow the user to add new
   %objectives/constraints
   cFormulationMenu = uimenu('Label','Update Formulation');
       uimNewObjectives = uimenu(cFormulationMenu,'Label','Add new objective(s)','enable',Boolean2Enabled(blValidLinProgramData),'Callback',{@AddObjectives,hWindReturn});
       uimNewConstraints = uimenu(cFormulationMenu,'Label','Add new constraint(s)','enable',Boolean2Enabled(blValidLinProgramData),'Callback',{@AddConstraints,hWindReturn});
    %define a local fontsize for the controls
    fontsizecntls = min([FontSize-4 12]);

    % Add the frames and controls
    %Frames

    %Frame to contain the tabs         
    ControlFrame = uipanel('Title','','FontSize',fontsizecntls,...
                   'Position',[xLeft+xWidth+lMarginPanel 0.005 xPanelWidthNorm 1-0.005],'BackgroundColor',get(hWindReturn,'Color'),'BorderType','none'); %[xLeft+xWidth+lMarginPanel 0.025 xPanelWidthNorm yBottom+yHeight]

    %Tabs
    hTabs = zeros(3,1);
    hTabButtons = zeros(3,1);
    sTabTips = {'' 'Tools to interactively filter and generate new alternatives' 'Tools to adjust the parallel coordinate plot display'};

    for i=1:3
       hTabs(i) = uipanel('Title','','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',ControlFrame,...
                 'Position',[0.02 0.025 .98 1-0.08]);
       hTabButtons(i) = uicontrol('Parent',ControlFrame,'Style', 'pushbutton','ToolTipString',sTabTips{i},'units','normalized','String', ButtonText{i},...
            'Position', [.05+.31*(i-1) .92 .31 .035],'fontsize', fontsizecntls);
       if i~=lStartTab
           set(hTabs(i),'visible','off');
       else
           set(hTabButtons(i),'BackgroundColor','white');
       end
    end    

    %set the callback functions for the tab buttons
    for i=1:3
        set(hTabButtons(i),'Callback', {@ShowTab,i,hTabs,hTabButtons,get(hTabButtons(2),'BackgroundColor')});
    end

    %Get the position of the tab in pixels
    %get(hTabs(1),'Position')
    set(hTabs(1),'Units','pixel');
    vPosTab = get(hTabs(1),'Position');
    set(hTabs(1),'Units','normalized');
    %get(hTabs(1),'Position')
     
    %The vertical position above the bottom of the tab box to place the
    %Group Decisions Checkbox
    lTopLine = vPosTab(4) - 55;
    lNearTop = vPosTab(4)-220; %500;    
    lTopColor = lTopLine; %Vertical position above the  bottom of the tab  box to place the Ramp Color checkbox
    lNearTopInteract = lTopLine; %570+70;
    lTopSlider = 241; %vPosTab(4) - 270; %565; %615;
    
    set(hTabs(1),'Units','Characters');
    hTabPos = get(hTabs(1),'Position');
    set(hTabs(1),'Units','normalized');
    lGSPanHeight = 11; %characters
    lGSPanPos = [1 hTabPos(4)-lGSPanHeight-9 hTabPos(3)-2 lGSPanHeight];
    
    if any(lGSPanPos<0)
       error('Figure height is too small. Either: 1) use a smaller font size, or 2) resize so Generate New Alterantives box on Interact tab is visible');
       lGSPanPos(lGSPanPos<0) = 0;
    end
    
    %Box and controls for Generate Solutions box
    GenSols = uipanel('Title','Generate New Alternatives','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',hTabs(1),...
                 'units','characters','Position',lGSPanPos); %[0.02 0.65 1-0.04 .2][7 vPosTab(4)-243 267 137]);
             % Bottom and height positions are in absolute characters ti specify absolute.
             % Set left and width position in normalized [0.02 0.70 1-0.04 .24](because the
             % parent Control Panel is in absolute characters).
        set(GenSols,'units','normalized');
        lGSPanPos = get(GenSols,'Position');
        lGSPanPos(1) = 0.02;
        lGSPanPos(3) = 1-0.04;
        set(GenSols,'Position',lGSPanPos);
        %Change to pixel for setting the controls within
        set(GenSols,'units','pixel');            
     
    %Box and controls for Filter Existing Alternatives
        %Position the Filter Existing Alternatives panel 0.02 normalized units below the
        %prior Move/Delete Panel. Height is the smaller of either the
        %remaining height of the parent Control Panel or box to tightly frame the group controls
        
    lFEAPanHeight = min([abs(lGSPanPos(2) - 0.02); lGSPanPos(4)*(n + 4)/7]);
    lFEAPanPos = [0.02 lGSPanPos(2) - 0.02 - lFEAPanHeight 1-0.04 lFEAPanHeight];
    
    if any(lFEAPanPos<0)
       error('Figure height is too small. Either: 1) use a smaller font size, or 2) resize so Filter Existing Alternatives box on Interact tab is visible');
       lFEAPanPos(lFEAPanPos<0) = 0;
    end
   
    FiltSols = uipanel('Title','Filter Existing Alternatives','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',hTabs(1),...
                 'units','normalized','Position',lFEAPanPos); %[0.02 0.01 1-0.04 .65-0.02]
             
     set(FiltSols,'Units','pixel');
     vFiltPos = get(FiltSols,'position');
     set(FiltSols,'Units','normalized');
   
    %Controls for the Interact Tab
    
    lblSliderValue = zeros(1,n);
    txtSliderValue = zeros(1,n);

    % convert steps to percentage of slider range
    %vStepUse = [(vMaxResult(:,1:nO)-vMinResult(:,1:nO))./100 vStep];
    vSliderStep = vStep./(vMaxResult - vMinResult);

    for i=1:n
        if (vSliderStep(i)>1) || isnan(vSliderStep(i))
            vSliderStep(i)=1;
        end
    end

    %vSliderStep

    %checkboxes at bottom of axes 
    cbChecks = zeros(1,n);
    
    %[vMinResult; vMaxResult];
    
    %Position Checkboxes
    for i=1:n
        cbChecks(i) = uicontrol('Parent',Figure2,'Style', 'checkbox', 'Units','normalized','Position',[xLeft+(YAxisMargins(1)+i-1)*xWidth/(n-1+sum(YAxisMargins)) yBottom 0.02 0.02],'fontsize', fontsizecntls,'Value',vFixed(i),...
                            'Callback',{@AxisChecked,cbChecks,i});%
        set(cbChecks(i),'units','characters');
        vPosCb = get(cbChecks(i),'Position');
 
        %vertical offset down 1 character, width and height a standard
        %characters as well
        vPosCbNew = [vPosCb(1) - 3.2/2 vPosCb(2)-1.6  3.2 chHeight];
        set(cbChecks(i),'Position',vPosCbNew);
        set(cbChecks(i),'units','normalized');
    end
      
    sCount = 0;
    maxSliders = 15;
            
    %Check box for hide/show slider control box
    %cbShowSliderFrame = uicontrol('Parent',hTabs(2),'Style', 'checkbox', 'String', 'Show slider panel','Position',[10 275 150 20],'fontsize', fontsizecntls,'Callback',{@ShowSliderFrame,hTabs(2)},'Value',1);    
    
    sSetBeyond = sprintf('When sliders are visible, check to set a value\non an axis beyond the range shown by the slider');
    
    %cbHideSliders = uicontrol('Parent',FiltSols,'Style', 'checkbox', 'Value',1,'String', 'Hide sliders','Position',[8 vFiltPos(4)-45 106 20],'fontsize', fontsizecntls,'value',HideSliders); %,'Callback',{@HideSliders,sSlider}); [20 lTopSlider 165 20]
      
    [sSlider,vShowSliders] = RenderSliders(vMinPlot,vMaxPlot,vFixedValsPlotUnits',vSliderStep);
    
    NearOptTolerance = zeros(nO,1); %Handles to the near optimal tolerance parameter text inputs for each objective
    sNearOpt = sprintf('Tolerance is the fraction of the optimal objective function value and\n adds a constraint to the underlying optimization problem\nto require alternatives with objective functon values within\nthe specified tolerance. Enter a number > 1 for minimization problems and < 1 for maximization problems.');
    
    %Label and control to set the top set value axis
    lblTopAxis = uicontrol('Parent',FiltSols,'Style', 'text','String','1st Decision axis to list (#):', 'Position', [10 vFiltPos(4)-55-16*(sCount+1.75) 175 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    txtTopAxis = uicontrol('Parent',FiltSols,'Style', 'edit', 'Position', [175 vFiltPos(4)-55-16*(sCount+1.75)-3 30 20],'fontsize', fontsizecntls-2);
    set(txtTopAxis,'string',1);
    
    %Header labels
    lblSetVals = uicontrol('Parent',FiltSols,'Style', 'text','String','Set value', 'Position', [10 vFiltPos(4)-55-16*(sCount+3) 55 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    lblAxis = uicontrol('Parent',FiltSols,'Style', 'text','String','Axis name', 'Position', [70 vFiltPos(4)-55-16*(sCount+3) 110 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    if nO>1
        lblNearOptTol = uicontrol('Parent',FiltSols,'Style', 'text','String','Near opt. tol.', 'Position', [70+115 vFiltPos(4)-55-16*(sCount+3) 85 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    end
    
    for i=1:n
        sCount = sCount+1;
%        if vShowSliders(i) > 0 
            sVis = 'on';
            if vFiltPos(4)-55-17*(sCount+3) < 0
                sVis = 'off'; %Hide sliders that fall outside the box
            end
            
            txtSliderValue(i) = uicontrol('Parent',FiltSols,'Style', 'edit', 'Position', [10 vFiltPos(4)-55-17*(sCount+3) 55 15],'fontsize', fontsizecntls-2, ...
                    'Visible',sVis,'Callback', {@SetSliderValue,sSlider,mTransformToOrig,i});

            lblSliderValue(i) = uicontrol('Parent',FiltSols,'Style', 'text','String',vAxisLabelsAll{i}, 'Position', [70 vFiltPos(4)-55-17*(sCount+3) 110 15], ...
                'Visible',sVis,'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
            %set(txtSliderValue(i),'String',sprintf('%.2f',get(sSlider(i),'Value')));
            set(txtSliderValue(i),'String',sprintf('%.2f',ConvertPlot2DataUnits(get(sSlider(i),'Value'),mTransformToOrig(:,i))));
            if (nO>1) && (i<=nO)
               %Add a textbox for the near-optimal alternative here
               NearOptTolerance(i) = uicontrol('Parent',FiltSols,'Style', 'edit','ToolTipString',sNearOpt,'Visible',sVis, 'Position', [70+115 vFiltPos(4)-55-17*(sCount+1.75) 45 15],'fontsize', fontsizecntls-2);
               set(NearOptTolerance(i),'string',Tolerance(i));
            end
            
%        end
    end     
    %circle back and now set the callback function for the sliders
    for i=1:n
       set(sSlider(i),'Callback',{@SetTxtValues,txtSliderValue,mTransformToOrig,i})
       set(txtSliderValue(i),'Callback', {@SetSliderValue,sSlider,mTransformToOrig,i});
    end  
    set(txtTopAxis,'Callback', {@SetAxesInFilterAltsBox,txtSliderValue,lblSliderValue,nO,vFiltPos(4)-55-17*(3),17});
    
    %Controls to set all sliders to a particular solution
    % Set to specified record #
    lblSetAll = uicontrol('Parent',FiltSols,'Style', 'text','String','Set to:', 'Position', [10 vFiltPos(4)-45 44 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');%[10 lTopSlider-63 60 15]
    txtCurrRec = uicontrol('Parent',FiltSols,'Style', 'edit','Position', [55 vFiltPos(4)-45 45 15],'fontsize', fontsizecntls-2,...
              'Callback', {@UpdateSetToControls,hWindReturn,0,0});%[10 lTopSlider-63 60 15]
    set(txtCurrRec,'String',lCurrRecord);
    
    sTotCount = sprintf('of %d',m);
    lblOfX = uicontrol('Parent',FiltSols,'Style', 'text','String',sTotCount, 'Position', [105 vFiltPos(4)-45 65 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');%[10 lTopSlider-63 60 15]

    %Buttons to advance records
    recFarLeft = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to the first record','String', '<<',...
            'Position', [55 vFiltPos(4)-65 20 15],'fontsize', fontsizecntls,....
            'Callback', {@UpdateSetToControls,hWindReturn,-1,2}); 
    recLeft = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to prior record','String', '<',...
            'Position', [75 vFiltPos(4)-65 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,-1,1});
    recRight = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to next record','String', '>',...
            'Position', [95 vFiltPos(4)-65 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,1,1});
    recFarRight = uicontrol('Parent',FiltSols,'Style', 'pushbutton','ToolTipString','... to last record', 'String', '>>',...
            'Position', [115 vFiltPos(4)-65 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,1,2});
    
    if lOptGroupCurr == 0
        sEnable = 'off';
    else
        sEnable = 'on';
    end
    
    GroupToHighlightButton = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'String', sprintf('Highlight Group'),...
            'Position',[167 vFiltPos(4)-65 87 40] ,'fontsize', fontsizecntls,'enable',sEnable,...
            'Callback', {@SetSliderValsToRecord,lHighRecord},'ToolTipString','Set record number and values to the highlighted group');
            % 'Callback', {@SetAllSliderVals,sSlider,txtSliderValue,mTransformToOrig,2,mPlot(OptGroupInds,:)}); %[150 lTopSlider-65 70 20]

    %Tool tips
    sAllowDeviation = sprintf('Is the fraction of a graph tick mark and used to define\n when an alternative is within the allowable deviation of a\nslider setting for a checked axis');
    
    %Label and Textbox for user to enter near-optimal tolerance if there is
    %only one
    NOToleranceLabel = uicontrol('Parent',hTabs(1),'Style', 'text','String','Near Optimal Tolerance (fraction of optimal):', 'Position', [10 lNearTopInteract-20 195 40],'fontsize', fontsizecntls,'BackgroundColor','white');
    NearOptToleranceAll = uicontrol('Parent',hTabs(1),'Style', 'edit','ToolTipString',sNearOpt, 'Position', [215 lNearTopInteract-20+10 45 20],'fontsize', fontsizecntls);
    
    if nO==1
        NearOptTolerance(1) = NearOptToleranceAll;
        set(NearOptTolerance,'String',Tolerance(1));
    else
        %Multi-objective, the NearOptToleranceAll entrywill update all the text
        %boxes
        
    end
    
    lblAllowableDeviation = uicontrol('Parent',hTabs(1),'Style', 'text','String','Allowable deviation (fraction of tick):', 'Position', [45 lNearTopInteract-55 145 30],'fontsize', fontsizecntls-2,'BackgroundColor','white');
    txtAllowableDeviation = uicontrol('Parent',hTabs(1),'Style', 'edit','ToolTipString',sAllowDeviation, 'Position', [215 lNearTopInteract-50 45 20],'fontsize', fontsizecntls-2,'Callback',@TriggerUpdateHighlightTraces);
    set(txtAllowableDeviation,'String',AllowableDeviation);
    

    
    %Drop-down box for Generate Type
    lTopGen = 90;
    GenerateType = 1;
    GenerateTypeLabel = uicontrol('Parent',GenSols,'Style', 'text','String','Type:','HorizontalAlignment','left', 'Position', [10 lTopGen 40 20],'fontsize', fontsizecntls-2,'BackgroundColor','white');
    cbGenerateType = uicontrol('Parent',GenSols,'Style','popup','String','One alternative|Maximum extents|Random sample|Enumerate all (MIPs)','Position',[59 lTopGen 125 25],'fontsize',fontsizecntls-2,'Value',GenerateType,'Callback',{@SetFromGenBoxes});
    set(cbGenerateType,'value',GenType);
    
    GenerateMethodLabel = uicontrol('Parent',GenSols,'Style', 'text','String','Engine:','HorizontalAlignment','left','Position', [10 lTopGen-30 45 20],'fontsize', fontsizecntls-2,'BackgroundColor','white');   
       
    sUsingTip = sprintf('Data - Query data already on the plot\nMatlab - Run Matlab linprog function using constraint system defined by\n          the ProbForm structure and (possibly) cFunc\nGAMS - Run the specified GAMS file'); 
    cbGenerateMethod = uicontrol('Parent',GenSols,'Style','popup','ToolTipString',sUsingTip,'String','Data|MATLAB (LP matrix)|GAMS','Position',[59 lTopGen-30 125 25],'fontsize',fontsizecntls-2,'Value',GenMethod,'Callback',{@SetFromGenBoxes});
 
    
    lblNumSamples = uicontrol('Parent',GenSols,'Style', 'text','String','# Samples:', 'Position', [10 lTopGen-55 75 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    txtNumSamples = uicontrol('Parent',GenSols,'Style', 'edit', 'Position', [80 lTopGen-55 45 20],'fontsize', fontsizecntls-2);
    lblObjSamples = uicontrol('Parent',GenSols,'Style', 'text','String','Obj. Samples (%):', 'Position', [130 lTopGen-55 120 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    sObjSamples = 'Percent of samples to stratify along objective function axes';
    txtObjSamples = uicontrol('Parent',GenSols,'Style', 'edit', 'Position', [240 lTopGen-55 30 20],'ToolTipString',sObjSamples,'fontsize', fontsizecntls-2);

    set(txtNumSamples,'String',num2str(NumSamples));
    set(txtObjSamples,'String',num2str(ObjSamplesPercent));
    %ResampleButton = uicontrol('Parent',hTabs(2),'Style', 'pushbutton', 'String', 'Resample',...
    %        'Position', [140 lTopSlider-30 90 20],'fontsize', fontsizecntls,'visible','off');
 
    lblGamsFile = uicontrol('Parent',GenSols,'Style', 'text','String','GAMS File:', 'Position', [10 lTopGen-80 75 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    sGamsFileTip = 'Full path and file name to GAMS file';
    txtGamsFile = uicontrol('Parent',GenSols,'Style', 'edit', 'ToolTipString',sGamsFileTip,'Position', [80 lTopGen-80 170 20],'fontsize', fontsizecntls-2);
    set(txtGamsFile,'String',sGamsFile);
    %RunGamsButton = uicontrol('Parent',hTabs(2),'Style', 'pushbutton', 'String', 'Run GAMS',...
    %        'Position', [180 lTopSlider-60 90 20],'fontsize', fontsizecntls,'visible','off');

    sGenSolsTip = sprintf('Generate new alternatives from the specifiednear-optimal tolerance\nand existing filtered alternatives');   
    GenerateButton = uicontrol('Parent',GenSols,'Style', 'pushbutton', 'String', 'Generate','ToolTipString',sGenSolsTip,...
            'Position', [189 lTopGen-25 75 40],'fontsize', fontsizecntls);
 
        
    
    %Check box for show/hide all controls (print view) -- goes above main frame
    %cbShowControls = uicontrol('Parent',ControlToggleFrame,'Style','checkbox', 'String', 'Show all controls','Position',[5 7 145 20],'fontsize', fontsizecntls,'Callback',{@ShowAllControls,ControlFrame,[plot2 ax2],cbChecks,sSlider,hTabContainer,xPanelWidthNorm},'Value',1);

    %Controls for the File Tab
    if 0
    %Button to load pareto optimal solutions
    ParetoCheck = uicontrol('Parent',hTabs(1),'Style', 'checkbox', 'String', 'Add','Position',[100 lTopLine-88 60 20],'fontsize', fontsizecntls);
    ParetoButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load Pareto',...
            'Position', [10 lTopLine-88 90 20],'fontsize', fontsizecntls,...
            'Callback', {@LoadPareto,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep,sSlider,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,FontSize,iOpt,mActCat,ParetoCheck});        % Pushbutton string callback

    %Button to permuate subregion
    TestButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Permute Subregion',...
            'Position', [28 lTopLine-110 128 20],'fontsize', fontsizecntls,...
            'Callback', {@TestAllRegion,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep,sSlider,mTransformToOrig,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,FontSize,iOpt,mActCat,vGroup,1});
    LoadTestResults = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load Results',...
            'Position', [28 lTopLine-135 128 20],'fontsize', fontsizecntls,...
            'Callback', {@TestAllRegion,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep, sSlider,mTransformToOrig,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,FontSize,iOpt,mActCat,vGroup,0});     
        
        
    %Button to plot the subregion
    PlotButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Plot Outer',...
            'Position', [28 lTopLine-160 128 20],'fontsize', fontsizecntls,...
            'Callback', {@RunButton,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,1,vararg_curr{:}});        % Pushbutton string callback
                                       % that calls a MATLAB function                                        
                                       
    %Button to reload gams output
    LoadGAMSButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load GAMS Output',...
            'Position', [28 lTopLine-185 128 20],'fontsize', fontsizecntls,...
            'Callback', {@RunButton,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,0,vararg_curr{:}});        % Pushbutton string callback
        %'Callback', {@RunButton,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep, sSlider,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,FontSize,iOpt,mActCat,vGroup,0});        % Pushbutton string callback
                                       % that calls a MATLAB function
    end
                                       
    %Controls for the Color Ramp Tab    

    %Controls for ramping color on traces along one axis
    CurrColorRamp = 0;
    CurrDirection = 1;
    sRampTip = sprintf('Ramp line colors on first checked axis\nin specified direction over specified # of classes\n(e.g., light to dark)');
    sDirectionTip = sprintf('Direction of color ramp:\nAscend: light to dark (small to large values)\nDescend: dark to light (small to large values)');
    cbRampColor = uicontrol('Parent',hTabs(2),'Style', 'checkbox', 'String', 'Ramp color','ToolTipString',sRampTip,'Position',[15 lTopColor 250 20],'fontsize', fontsizecntls,'Callback',{@RampColor},'Value',CurrColorRamp);
    lblNumClasses = uicontrol('Parent',hTabs(2),'Style', 'text', 'String','# Color classes:','Position',[15 lTopColor-1*25 130 20],'fontsize', fontsizecntls,'BackgroundColor','white','HorizontalAlignment','Left');
    txtNumClasses = uicontrol('Parent',hTabs(2),'Style', 'edit', 'Position',[135 lTopColor-1*25 30 20],'fontsize', fontsizecntls,'String','10');
    
    lblDirection = uicontrol('Parent',hTabs(2),'Style','text','String','Direction:','Position',[15 lTopColor-2*25 130 20],'fontsize', fontsizecntls,'BackgroundColor','White');
    cbDirection = uicontrol('Parent',hTabs(2),'Style','popup','ToolTipString',sDirectionTip,'String','Ascend|Descend','Position',[135 lTopColor-2.5*20 75 20],'fontsize', fontsizecntls-2,'Callback',{@RampColor},'Value',CurrDirection);
    
    lNumRowsForColorChecks = floor((lTopColor-45)/20);
    
    
    %Controls for the Display Tab
        set(hTabs(3),'Units','Characters');
        hTabPos = get(hTabs(3),'Position');
        set(hTabs(3),'Units','normalized');
        
        %Formating panel
        lFormPanHeight = 6; %characters
        FormatPanel = uipanel('Title','Formatting','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',hTabs(3),...
                 'units','characters','Position',[1 hTabPos(4)-lFormPanHeight-3 hTabPos(3)-2 lFormPanHeight]);
             % Bottom and height positions are in absolute characters ti specify absolute.
             % Set left and width position in normalized [0.02 0.70 1-0.04 .24](because the
             % parent Control Panel is in absolute characters).
             % 
        set(FormatPanel,'units','normalized');
        lFormPanPos = get(FormatPanel,'Position');
        lFormPanPos(1) = 0.02;
        lFormPanPos(3) = 1-0.04;
        set(FormatPanel,'Position',lFormPanPos);
               
        %Change to pixel for setting the controls within
        set(FormatPanel,'units','pixel');
        lFormPanPosPx = get(FormatPanel,'Position');
    
        lTopFormat = lFormPanPosPx(4) - 15;
        lblFontSize = uicontrol('Parent',FormatPanel,'Style', 'text','String','Font Size:', 'Position', [10 lTopFormat-35 75 20],'HorizontalAlignment','Left','fontsize', fontsizecntls,'BackgroundColor','white');
        txtFontSize = uicontrol('Parent',FormatPanel,'Style', 'edit', 'Position', [90 lTopFormat-35 45 20],'fontsize', fontsizecntls);
        set(txtFontSize,'String',num2str(FontSize));
        
        %Move/Delete Axes Panel
        lMAPanHeight = min([hTabPos(4)-lFormPanHeight-3; 14]);   %characters
        
        MoveAxesPanel = uipanel('Title','Move/Delete Axes','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',hTabs(3),...
                 'units','characters','Position',[1 hTabPos(4)-lFormPanHeight-lMAPanHeight-4 hTabPos(3)-2 lMAPanHeight]);
             % Bottom and height positions are in absolute characters ti specify absolute.
             % Set left and width position in normalized [0.02 0.70 1-0.04 .24](because the
             % parent Control Panel is in absolute characters).
             % 
        set(MoveAxesPanel,'units','normalized');
        lMAPanPos = get(MoveAxesPanel,'Position');
        lMAPanPos(1) = 0.02;
        lMAPanPos(3) = 1-0.04;
        set(MoveAxesPanel,'Position',lMAPanPos);
               
        %Change to pixel for setting the controls within
        set(MoveAxesPanel,'units','pixel');
    
        lTopMoveAxis = 250;
        
        
    %Check box for show decision variable labels
    sShowGrouping = sprintf('Show a table below the axes labels\ncomprised of inputs from mGroupData\n that separates & colors decision axes\nby major and minor categories');
    sHideChecks = sprintf('Hide the check boxes below each axis');

    %Label and Textbox for move first checked axis left or right
    MoveAxisLabel = uicontrol('Parent',MoveAxesPanel,'Style', 'text','String','Move checked axes:', 'HorizontalAlignment','center','Position', [20 lTopMoveAxis-6.75*20 100 40],'fontsize', fontsizecntls,'BackgroundColor','white');

    sFarLeft = ['... to the far left'];
    sFarRight = ['... to the far right'];
    sLeft = ['... one position to the left'];
    sRight = ['... one position to the right'];
    
   
    
    FarLeftButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton', 'ToolTipString',sFarLeft,'String', '<<',...
            'Position', [125 lTopMoveAxis-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults, nO, 1, -1,NearOptTolerance,txtAllowableDeviation,cbChecks, sSlider,mTransformToOrig,vararg_out{:}});
        %     'Callback', {@MoveAxis,hWindReturn,mResults, nO, 1, -1,NearOptTolerance, cbChecks, vFixedVals, vStep, sSlider,vObjLabels,vXLabels, vXLabelsShort, yAxisLabels, FontSize, iOpt, mActCat, vGroup});
 
    LeftButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton', 'ToolTipString',sLeft,'String', '<',...
            'Position', [145 lTopMoveAxis-6.75*20+15 20 20],'fontsize', fontsizecntls);
        %    'Callback', {@MoveAxis,mResults,nO, 0, -1,NearOptTolerance,txtAllowableDeviation, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});
    RightButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton', 'ToolTipString',sRight,'String', '>',...
            'Position', [165 lTopMoveAxis-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults,nO, 0, 1,NearOptTolerance,txtAllowableDeviation, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});
    FarRightButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton','ToolTipString',sFarRight, 'String', '>>',...
            'Position', [185 lTopMoveAxis-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults,nO, 1, 1,NearOptTolerance,txtAllowableDeviation, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});

         
   %checkbox that determines how axes are re-ordered
    sReorder = sprintf('Reorder axes on the plot\nby the selected feature.');
    sByCat = '... by category specified in the input mActCat';
    sByFin = sprintf('... by values on the axis\n(e.g., first plot axes with singular (positive) values,\nnext axes with positive ranges,\n last axes that all are zero)');
    cbByCat = uicontrol('Parent',MoveAxesPanel,'Style', 'checkbox','ToolTipString',sByCat, 'String', 'By Category','Position',[90 lTopMoveAxis-8.5*20+5 145 20],'fontsize', fontsizecntls);
    cbReorder = uicontrol('Parent',MoveAxesPanel,'Style', 'checkbox','ToolTipString',sByFin, 'String', 'Defin. actions 1st','Position',[90 lTopMoveAxis-9.5*20+5 145 20],'fontsize', fontsizecntls);%    


    %Button to reorder axes
    ReorderButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton','ToolTipString',sReorder,...
            'Position', [15 lTopMoveAxis-9.5*20 70 50],'fontsize', fontsizecntls,...
            'Callback', {@Reorder,mResults,nO,hWindReturn});
                    %'Callback', {@Reorder,mResults, nO,NearOptTolerance, cbChecks, sSlider,mTransformToOrig,cbReorder,cbByCat,vararg_out{:}});        % Pushbutton string callback

    %Wrap the text in ReorderButton
    cText = {'Reorder axes'};
    [outstring,newpos] = textwrap(ReorderButton,cText);
    set(ReorderButton,'String',outstring); %,'Position',newpos);        
                 

        
    %Button to reorder decision variables
    PruneButton = uicontrol('Parent',MoveAxesPanel,'Style', 'pushbutton','ToolTipString','Remove selected axes', ...  %'HorizontalAlignment','center',
            'Position', [15 lTopMoveAxis-12*20 70 50],'fontsize', fontsizecntls); %set the callback function after define AllowableDeviation
        
    %Wrap the text in PruneButton
    cText = {'Delete axes'};
    [outstring,newpos] = textwrap(PruneButton,cText);
    set(PruneButton,'String',outstring); %,'Position',newpos);    
        
    %checkboxes and buttons to prune axes
    cbPruneChecked = uicontrol('Parent',MoveAxesPanel,'Style', 'checkbox', 'String', 'Checked','Position',[90 lTopMoveAxis-11*20+5 145 20],'fontsize', fontsizecntls);
    cbPruneZeros = uicontrol('Parent',MoveAxesPanel,'Style', 'checkbox', 'String', 'All zero','Position',[90 lTopMoveAxis-12*20+5 145 20],'fontsize', fontsizecntls);%    

    GroupTracesPanel = 0;

    %Trace Group checkboxes to toggle their visibility
    if nU>0
        
        %Position the Group Traces Panel 0.02 normalized units below the
        %prior Move/Delete Panel. Height is the smaller of either the
        %remaining height of the parent Control Panel or box to tightly frame the group controls
        
        lGTPanHeight = min([lMAPanPos(2) - 0.02; lMAPanPos(4)*(nU + 6)/6]);
               
        GroupTracesPanel = uipanel('Title','Group Traces','FontSize',fontsizecntls,...
                 'BackgroundColor','white','parent',hTabs(3),...
                 'units','normalized','Position',[0.02 lMAPanPos(2) - 0.01 - lGTPanHeight 1-0.04 lGTPanHeight]);
             
        %Set units to pixel to position controls within
        set(GroupTracesPanel,'units','pixel');
        gtPos = get(GroupTracesPanel,'position');
        lTopGroupTraces = gtPos(4)-35;           
             
        %Build the popup selection options for groups
        sGroupOptions = ['-- none --|',strjoin(vUnique','|')];
        lblHighGroup = uicontrol('Parent',GroupTracesPanel,'Style', 'text', 'String','Highlight Group:','Position',[15 lTopGroupTraces-(1)*20 120 20],'fontsize', fontsizecntls,'ForegroundColor',[0 0 0]);
        
        cbHighlightGroup = uicontrol('Parent',GroupTracesPanel,'Style', 'popup', 'String', sGroupOptions,'Position',[140 lTopGroupTraces-(1)*15 110 20],'fontsize', fontsizecntls-2,'ForegroundColor',[0 0 0],'value',lOptGroupCurr+1, ...
                'Callback', {@HighlightGroup,mResults,nO,hWindReturn});         

        
        cbGroupChecks = zeros(nU,1);
        txtGroupOrders = zeros(nU,1);
        txtGroupNames = zeros(nU,1);
        txtGroupThicks = zeros(nU,1);

        sReorderGrp = sprintf('Reorder groups by number in the Ord column');
        sMergeGrp = 'Combine checked groups into one group (first checked group)';
        
        ShowGroupsBut = uicontrol('Parent',GroupTracesPanel,'Style', 'pushbutton','String','Show Checked Groups', 'HorizontalAlignment','center', 'Position', [15 lTopGroupTraces-2.5*20+5 195 20],'fontsize', fontsizecntls, ...
              'Callback', {@ShowGroups,hWindReturn,2});
        RemoveGroupsBut = uicontrol('Parent',GroupTracesPanel,'Style', 'pushbutton','String','Delete Checked Groups', 'HorizontalAlignment','center', 'Position', [15 lTopGroupTraces-3.5*20+5 195 20],'fontsize', fontsizecntls, ...
              'Callback', {@RemoveGroups,mResults,nO,hWindReturn,2});
        ReorderGroupsBut = uicontrol('Parent',GroupTracesPanel,'Style', 'pushbutton','ToolTipString',sReorderGrp,'String','Reorder Groups', ....
                  'HorizontalAlignment','center', 'Position', [15 lTopGroupTraces-4.5*20+5 195 20],'fontsize', fontsizecntls);
        MergeGroupsBut = uicontrol('Parent',GroupTracesPanel,'Style', 'pushbutton','ToolTipString',sMergeGrp,'String','Merge Checked Groups', ....
                  'HorizontalAlignment','center', 'Position', [15 lTopGroupTraces-5.5*20+5 195 20],'fontsize', fontsizecntls);

              
        %Row headers for the group entries
        lblGroupOrder = uicontrol('Parent',GroupTracesPanel,'Style','text','String','Ord.','HorizontalAlignment','left','Position',[35 lTopGroupTraces-7*20+5 28 20],'fontsize', fontsizecntls);
        lblGroupName = uicontrol('Parent',GroupTracesPanel,'Style','text','String','Name','HorizontalAlignment','left','Position',[65 lTopGroupTraces-7*20+5 145 20],'fontsize', fontsizecntls);
        lblGroupThickness = uicontrol('Parent',GroupTracesPanel,'Style','text','String','Thick.','HorizontalAlignment','left','Position',[212 lTopGroupTraces-7*20+5 45 20],'fontsize', fontsizecntls);
             
              
        for i=1:nU
            if i==lOptGroupCurr
                cColor = mHighlightColor;
            else
                lIndToUse = 1+mod(i-1,size(mColors,1));
                cColor = mColors(lIndToUse,1,3,:);
            end

            cbGroupChecks(i) =  uicontrol('Parent',GroupTracesPanel,'Style', 'checkbox', 'String', '','Position',[15 lTopGroupTraces-(7+i)*20+5 20 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'value',vShowGroup(i));
            txtGroupOrders(i) =  uicontrol('Parent',GroupTracesPanel,'Style', 'edit','ToolTipString','Number to specify group order in plot','String', num2str(i),'Position',[35 lTopGroupTraces-(7+i)*20+5 28 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'BackgroundColor','white');
            txtGroupNames(i) =  uicontrol('Parent',GroupTracesPanel,'Style', 'edit','ToolTipString','Label for group','String', vUnique(i),'Position',[65 lTopGroupTraces-(7+i)*20+5 145 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'BackgroundColor','white','HorizontalAlignment','left', ...
                                    'Callback',{@RenameGroup,i});
                                
            txtGroupThicks(i) =  uicontrol('Parent',GroupTracesPanel,'Style', 'edit','ToolTipString','Line thickness for group','String', vGroupThick(i),'Position',[212 lTopGroupTraces-(7+i)*20+5 45 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'BackgroundColor','white', ...
                                    'Callback', {@ShowGroups,hWindReturn,2});

        end

        cbGroupsOneColor = uicontrol('Parent',GroupTracesPanel,'Style', 'checkbox','ToolTipString','Plot all groups with a single color','String', 'All Groups Same Color','Position',[15 lTopGroupTraces-(8+nU)*20+5 195 20],'fontsize', fontsizecntls,'ForegroundColor',[0 0 0],'value',0);


    else
        cbGroupChecks=[];
        txtGroupOrders=[];
        txtGroupNames=[];
        txtGroupThicks=[];
    end
 
    %Build a cell array of handles to the controls we'll need to query in
    %callback functions. Organize by single-value text Box inputs, Axis
    %inputs, Group inputs, and other UI settings.
    control_handles = {'Tolerance' NearOptTolerance 'AllowableDeviation' txtAllowableDeviation 'NumSamples' txtNumSamples ...
                        'ObjSamples' txtObjSamples ...
                        'Tabs' hTabs 'AxisChecked' cbChecks 'Sliders' sSlider 'txtSliderValue' txtSliderValue ...
                        'GroupChecks' cbGroupChecks 'GroupOrders' txtGroupOrders 'GroupNames' txtGroupNames 'GroupThicks' txtGroupThicks  ...
                        'HideSliders' uimHideSliders 'HideCheckboxes' uimHideChecks 'ByCat' cbByCat 'Reorder' cbReorder 'PruneChecked' cbPruneChecked ...
                        'PruneZeros' cbPruneZeros 'AllGroupsSameColor' cbGroupsOneColor 'HighlightGroup' cbHighlightGroup 'GamsFile' txtGamsFile ...
                        'ColorRamp' cbRampColor 'RampDirection' cbDirection 'NumClasses' txtNumClasses 'cbGenType' cbGenerateType 'cbGenMethod' cbGenerateMethod ...
                        'CurrentRecord' txtCurrRec 'TotalRecords' lblOfX 'GroupToHightlightBut' GroupToHighlightButton 'HideControls' uimShowAll ...
                        'txtFontSize' txtFontSize 'ShowCurrRecord' uimShowCurrRec 'ShowFilteredAlts' uimShowFiltered, 'ShowObjsDiffColor' uimShowObjsDiffColor...
                        'ShowGroupLabels' uimGroupLabels,'ShowInsetPlot',uimShowInset,'LegendOptions',uimLegendOptions}; 

    %Also set as an app variable
    setappdata(hWindReturn,'hControls',control_handles);                
                    
    set(FarLeftButton,'Callback', {@RearrangeAxes,1,-1,hWindReturn,mResults,nO});
    set(LeftButton,'Callback', {@RearrangeAxes,0,-1,hWindReturn,mResults,nO});
    set(RightButton,'Callback', {@RearrangeAxes,0,1,hWindReturn,mResults,nO});
    set(FarRightButton,'Callback', {@RearrangeAxes,1,1,hWindReturn,mResults,nO});
    
    set(PruneButton,'Callback', {@RearrangeAxes,0,0,hWindReturn,mResults,nO});
    %set(ResampleButton,'Callback',{@Resample,NearOptTolerance,txtAllowableDeviation,txtNumSamples,cbGroupChecks,cbChecks,sSlider,mTransformToOrig,mResults,nO,hWindReturn});
    %set(ResampleButton,'Callback',{@Resample,hWindReturn,control_handles,mTransformToOrig,mResults,nO});
    %set(RunGamsButton,'Callback',{@EnumerateWithGams,hWindReturn,control_handles,mTransformToOrig,mResults,nO});
    set(GenerateButton,'Callback',{@GenerateNewAlts,hWindReturn});
    
    set(RemoveGroupsBut,'Callback',{@RemoveGroups,mResults,nO,hWindReturn,2});
    set(MergeGroupsBut,'Callback',{@RemoveGroups,mResults,nO,hWindReturn,1});
    set(ReorderGroupsBut,'Callback',{@ReorderGroups,mResults,nO,hWindReturn})
    set(txtFontSize,'Callback',{@ReorderGroups,mResults,nO,hWindReturn})
    %New callback to allow calling nearoptplot2 at the end
    %set(LoadTestResults,'Callback',{@TestAllRegion,hWindReturn,mResults,nO,control_handles,mTransformToOrig,0});
    
    
   %Callbacks for menu items
   %set(uimLoadData,'Callback',{@RunButton,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,0,vararg_curr{:}});
   %set(uimLoadData,'Callback',{@TestAllRegion,hWindReturn,mResults,nO,0});
   set(uimShowAll,'Callback',{@ShowAllControls,ControlFrame,[plot2 ax2],cbChecks,sSlider,hTabContainer,xPanelWidthNorm+rMarginPanel});  
     set(uimHideSliders,'Callback',{@ToggleSliders}); %{@callHideSliders,sSlider});
     set(uimSetBeyond,'Callback',{@ToggleMenuItem});
     set(uimCheckAll,'Callback',{@CheckAllBoxes,cbChecks})
     set(uimHideChecks,'Callback',{@ToggleCheckBoxVisibility,cbChecks,hTexts,hTabContainer});
   set(uimGroupLabels,'Callback',{@GroupDecisionLabels,hTexts,hTabContainer,nO,plot2});
     set(uimShowCurrRec,'Callback',{@ToggleCurrMenuItem});
     set(uimShowFiltered,'Callback',{@ToggleCurrMenuItem});
     set(uimShowObjsDiffColor,'Callback',{@ToggleShowObjs});
     set(uimShowInset,'Callback',{@ToggleInset});
   for i=1:length(uimTab)
       set(uimTab(i),'Callback',{@ShowTab,i,hTabs,hTabButtons,get(hTabButtons(2),'BackgroundColor')});
   end


%Set Figure units to characters to allow multi-platform plotting
set(hWindReturn,'Units','normalized');
%Set all other object units to normalized to allow multi-platform plotting
set(cbChecks,'Units','normalized');
set(sSlider,'Units','normalized');
set(hTexts,'Units','normalized');
for i=1:3
    set(get(hTabs(i),'Children'),'Units','normalized');
end
hPanels = [GenSols FiltSols MoveAxesPanel GroupTracesPanel FormatPanel];
for i=1:length(hPanels)
    set(hPanels(i),'Units','normalized');
    set(get(hPanels(i),'Children'),'Units','normalized');
end

%% Data for the mouse event callback functions

%Create a dummy over- parallel coordinate plot to highlight mouse-overed
%solutions
% Find out which groups are checked
gDexes = GetAllChecked(cbGroupChecks);

hold on
MouseOverPCPDec=parallelcoords([],'color',mColors(gDexes(1),1,3,:)); % GreenToMagentaRamp(1));
SubSpacePCPDec=parallelcoords([],'color',mColors(gDexes(1),1,2,:)); %GreenToMagentaRamp(3));

%Create the ramps 
[ColorRampPCPs,hColorChecks] = CreateEmptyColorRampPCPs(txtNumClasses);

if (nO>0) %overplot the objective axes in a different color
    MouseOverPCPObj=parallelcoords([],'color',mColors(gDexes(1),2,3,:)); %GreenToMagentaRamp(16));
    SubSpacePCPObj=parallelcoords([],'color',mColors(gDexes(1),2,2,:)); %GreenToMagentaRamp(13));
end

textHdl = text('Color', 'black', 'VerticalAlign', 'Bottom');

%Trigger callbacks for menu item settings that were passed as input
%arguements to the plot
if lShowGroupLabels == 1 %Trigger the function to show the table
    GroupDecisionLabels(uimGroupLabels,0,hTexts,hTabContainer,nO,plot2)
else %Just update the table position
    PositionGroupLabelsTable(hTabContainer, nO, hTexts, plot2)
end

hLegend = 0; %Handle for future legend

if ShowControls == 0 %Trigger the function to hide the controls and possibly the legend 
    ShowAllControls(uimShowAll,0,ControlFrame,[plot2 ax2],cbChecks,sSlider,hTabContainer,xPanelWidthNorm+rMarginPanel);
else
    UpdateLegend;
end
if HideSliders==0
   UpdateHighlightedTraces;
elseif HideCheckboxes==1
   ToggleCheckBoxVisibility(uimHideChecks,0,cbChecks,hTexts,hTabContainer)
end


    
%Create the cartesian inset graph
if nO==2
    [hInsetPlot,hBackBox] = LaunchInset(hWindReturn);
    if ShowInsetPlot
        ToggleInset(uimShowInset,0);
    end
end


%% Callback functions for mouse events

%Set the mouse call backs for moving the mouse and clicking to fix a value
set(hWindReturn,'WindowButtonMotionFcn', @MouseHoverCallback);
set(hWindReturn,'WindowButtonDownFcn', @MouseButtonDownCallback);
set(hWindReturn,'ResizeFcn',@ResizeCallback);  

%% The Mouse Call back functions

    function ResizeCallback(src,evt)
        %Callback when the figure is resized. Grab the current data and
        %control settings and replot the figure (to update the positioning)      
        %ReorderGroups(0,0,mResults,nO,hWindReturn);
        PositionGroupLabelsTable(hTabContainer, nO, hTexts, plot2);
    end

   function MouseButtonDownCallback(src,evt)
       %Callback for the mouse click event
       %Fix the selected axis at the value clicked. Highlights  all
       %lines that pass through this point and other prior fixed points
         
        axesHdl = plot2;
        % Grab the x & y axes coordinate where the mouse is
        mousePoint = get(axesHdl, 'CurrentPoint');
        mouseX = mousePoint(1,1);
        mouseY = mousePoint(1,2);

        %Use the X coordinate to find out which axis the mouse is near
        cI = round(mouseX);
        
         yrange = range(get(axesHdl, 'Ylim'));

        %abort if the mouse is not near the axis or the axis is checked
        %(already fixed)
        if (abs(mouseX - cI) > 0.2)  || (cI>n) || (cI<1) % || (mouseY < yrange(1)) || (mouseY>yrange(2))
            return
        end
        
        sMaxValuePU = get(sSlider(cI),'Max');
        sMinValuePU = get(sSlider(cI),'Min');   
        
        if mouseY>sMaxValuePU 
           %outside feasible range. Set to max
           mouseY=sMaxValuePU;
           warning('Clicked outside feasible range. Setting to max value')
        end
        if mouseY<sMinValuePU 
           %outside feasible range. Set to min
           mouseY=sMinValuePU;
           warning('Clicked outside feasible range. Setting to min value')
        end     
        
        %Register the click
        %Set the checkbox
        set(cbChecks(cI),'Value',1);
        %Set the slider; should also trigger the text box to update
        set(sSlider(cI),'Value',mouseY);
        SetTxtValues(sSlider(cI),[],txtSliderValue,mTransformToOrig,cI);
        UpdateHighlightedTraces;
        %ConvertPlot2DataUnits(mouseY,mTransformToOrig(:,cI),0)
   end

   function MouseHoverCallback(src,evt)
       %Callback for Mouse Hover
       %Undo prior highlighting and highlight the lines that pass through the point where the mouse currently is
       %These lines are limited to the first checked group as well as fixed
       %values set on other axis
       
        if any(~ishandle(cbGroupChecks))
            %Object we need doesn't exist. Ignore callback
            return
        end
    
        axesHdl = plot2;
        sShowObjsDiffColor = get(uimShowObjsDiffColor,'Checked');
        %delete the existing mouse-overed parallelcoords plot and create a
        %new dummy one
        delete(MouseOverPCPDec);
        if (nO>0) %%Plotting objectives in a different color
            delete(MouseOverPCPObj);
        end
        hold on
        
        %Find the groups to work with
        
        try
            OrderedGroups = vUnique; %unique(vGroup);
            [gDexes] = GetAllChecked(cbGroupChecks);
        catch
            %Some kind of error with cbGroupChecks; quick
           return; 
        end
            
            
            vLineStyleToUse = vLineStyle{gDexes(1)};
        
        %Use that group to see the colors
        MouseOverPCPDec = parallelcoords([],'color',mColors(gDexes(1),1,3,:)); %GreenToMagentaRamp(1));
        if (nO>0) %%Plotting objectives in a different color
            MouseOverPCPObj = parallelcoords([],'color',mColors(gDexes(1),2,3,:)); %GreenToMagentaRamp(16));
        end
        
        set(textHdl, 'String', '')
        % Grab the x & y axes coordinate where the mouse is
        mousePoint = get(axesHdl, 'CurrentPoint');
        mouseX = mousePoint(1,1);
        mouseY = mousePoint(1,2);

        %Use the X coordinate to find out which axis the mouse is near
        cI = round(mouseX);
        
        yrange = (get(axesHdl, 'Ylim'));

        %abort if the mouse is not near the axis or the axis is checked
        %(already fixed)
        if (abs(mouseX - cI) > 0.2)  || (cI>n) || (cI<1) || (mouseY < yrange(1)) || (mouseY>yrange(2))
            return
        end

        %Read in the variable arguments
        %[vFixedVals vGroup] = aGetField(varargin,{'vFixedVals' 'vGroup'}); 
        % Read in the AllowableDeviationVal
        if ShowControls==1
            [AllowableDeviationVal, AllowableDeviationValPU] = GetCheckAllowableDeviation(txtAllowableDeviation,'MouseMove',TickSize);

            %Find the rows of mResults that are checked and fixed
            [cRows] = ReturnRows(vGroup,OrderedGroups(gDexes),AllowableDeviationValPU,mResults,cbChecks,vFixedVals,sSlider,mTransformToOrig);
            %ErrorTolPU = ConvertPlot2DataUnits(ErrorTolDU,mTransformToOrig(:,cI),1);

            %Search along the cI-th column of mResults to find values that within
            %the ErrorTol of the y-mouse position
            RowsToHighLight=cRows(abs(mPlot(cRows,cI)-mouseY) <= AllowableDeviationValPU);
        else
            RowsToHighLight=[];
        end
        
        xrange = range(get(axesHdl, 'Xlim'));
        yrange = (get(axesHdl, 'Ylim'));
        
        set(textHdl, 'String',ThousandSep(ConvertPlot2DataUnits(mouseY,mTransformToOrig(:,cI),0)));
        set(textHdl, 'Position', [cI + 0.01*xrange, mouseY + 0.01*yrange])
        
        if ~isempty(RowsToHighLight)
            MouseOverPCPDec = parallelcoords(mPlot(RowsToHighLight,:),'color',mColors(gDexes(1),1,3,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize);
            hold on
            if (nO>0) %%Overplot objectives in a different color
                MouseOverPCPObj = parallelcoords(mPlot(RowsToHighLight,1:nO+1),'color',mColors(gDexes(1),2,3,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize,'Visible',sShowObjsDiffColor);
            end
           
        end
   end

    function UpdateHighlightedTraces
        % Updates the highlighted traces that are overplotted in darker
        % colors
        % 
        % The Source of traces depends on whether sliders are hidden or
        % visible:
        %   - Hidden: The current record (if Controls Menu=>Show current record checked)
        %   - Visible: Sub space defined by fixed values, checked axis, Type of generation, and method used to generate   
        %              Also updates the sliders that show the maximal
        %              extents of each axis
        %

        sSliderUse = sSlider;
        
        %delete the existing highlighted traces and create a
        %new dummy one
        delete(SubSpacePCPDec);
        if nO>0
            delete(SubSpacePCPObj);
        end

        %Find the groups to work with
        OrderedGroups = vUnique; %unique(vGroup);
        [gDexes] = GetAllChecked(cbGroupChecks);
        vLineStyleToUse = vLineStyle{gDexes(1)};

        hold on

        %Read whether the hide sliders and show current record menu items are checked is selected
        %Query the relevant controls and structures
        [vParams,hControls] = ReadControls('UpdateHighlightTraces',hWindReturn);
        [sGamsFile,vFixed,vFixedVals,ToleranceValue,lCurrRecord,LimitSource,lHideSliderValue,lShowCurrRecord,lShowFiltered,lShowObjsDiffColor] = aGetField(vParams,{'sGamsFile' 'vFixed' 'vFixedVals' 'Tolerance' 'CurrentRecord' 'GenerateMethod' 'HideSliders' 'ShowCurrRecord' 'ShowFilteredAlts' 'ShowObjsDiffColor'});
       
        sObsVis = Boolean2Enabled(lShowObjsDiffColor);
        
        %Error check validity of limit settings. If not valid, default to a
        %LimitSource setting of 1 (from Data records)     
        if (LimitSource==2) && (isempty(AMat) || isempty(Brhs) || isempty(OptSolRow) || isempty(NearOptConstraint))
            warning('Linear system not defined (ProbForm). Continuing with plot data');
            set(cbGenerateMethod,'Value',1);
            LimitSource = 1;
        end
        if (LimitSource==3) && ~exist(sGamsFile,'file')
            warning('Gams file does not exist. Continuing with plot data');
            LimitSource = 1;
        end
       
        n = length(cbChecks);
        mResultsPlot = [];
        mResultsCurrRec = [];
        lWidth = 1.25;
 
        if (lShowCurrRecord==1)
            %Show the current record
            mResultsCurrRec = mPlot(lCurrRecord,:);
            lWidth = 1.75;
        end
        if (lHideSliderValue==0) && (LimitSource==1)
            %Sliders visible, draw from exists of existing data
            [mResCompact,mResultsData,cRows] = ReturnDataExtents(hWindReturn,mResults,nO,1);           
            mResultsPlot = mResultsData;
        elseif (lHideSliderValue==0) && (LimitSource==2)       
            %Sliders visible, query min and max values using linear constraints defined in
            %ProbForm (problem formulation). Use maximum extent functions
            [ProbNew,ProbOld,AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,AllowableDeviationValue] = UpdateLPMatrix(hWindReturn,nO,2);
            
            [mResultsValUse,mResCompact] = maxextentind(ProbNew,struct('vFixedVals',vFixedVals(nO+1:end)','vFixed',vFixed(nO+1:end),'UnfixCurrX',1)); %'Aeq',Aeq,'Beq',Beq,
            
            if any(isnan(mResultsValUse))
                warning('Can not calculate maximum extents for current settings. Stopping');
                mResultsPlot = [];
            end
                        
            mResCompact = mResCompact';
            [mDF,nDF] = size(mResultsValUse);
            
            if (nDF < nD) && (nDF + sum(vFixed(nO+1:end)) == nD)
                mResCompact = CompactToFull(mResCompact,vFixed(nO+1:end),vFixedVals(nO+1:end)');
                mResultsValUse = CompactToFull(mResultsValUse,vFixed(nO+1:end),vFixedVals(nO+1:end)');
            end
            
                       
            if (nO>0) && (~isempty(cFunc))
                %Find the range of objective function values
                [mResultsValObjs,mResCompObjs] = maxextentind(ProbNew,struct('objfunc',cFuncFree')); %Algorithm simplex
                if nDF < nD
                    [mResultsValObjs] = CompactToFull(mResultsValObjs,vFixed(nO+1:end),vFixedVals(nO+1:end)');
                end
                mResultsValUse = [mResultsValUse;mResultsValObjs];

                %Calculate the objective function values for exsiting solutions           
                mResultsValUse = [mResultsValUse*cFunc' mResultsValUse];
                %Update the Compact Version
                mObjComp = [];
                for i=1:nO
                    mObjComp = [mObjComp mResultsValUse(mDF+2*[(i-1):i],i)];
                end
                mResCompact = [mObjComp mResCompact];
            else
                warning('Need to define cFunc to calculate objective function values')
            end            
            
            sprintf('Current Allowable Axes Ranges based on Slider/Checkbox Settings')
            {'Axis' 'Fixed'  'Min' 'Fixed Value' 'Max'}
            [num2cell([1:nO 1:nD]') num2cell(vFixed') num2cell(mResCompact(1,:)') num2cell(vFixedVals) num2cell(mResCompact(2,:)')]
            
            for i=1:n
                mResultsPlot(:,i) = ConvertPlot2DataUnits(mResultsValUse(:,i),mTransformToOrig(:,i),1);
            end           
        elseif (lHideSliderValue==0) && (LimitSource==3)
            %Use GAMS to return the limits
            %Change the directory/folder to the one the GAMS File is in
            [sPath,sFileName,sFileExt] = fileparts(sGamsFile);
            if ~strcmpi(sPath,'')
                cd(sPath);
            end

            RunMode = 2;
            ColInds = [1:n];      

            %Call GAMS to identify the maximum extents            
            [vObjs, mResultsInt, mResultsVal, uelsOut, vReturnFlag, mGamsStats, NumSolvs] = EnumNEIntSolsGams4(sGamsFile,vFixed(nO+1:nO+nD),vFixedVals(nO+1:nO+nD),vXLabelsShort,ToleranceValue(1),2,RunMode);

            NumFeasibleSols = sum(vReturnFlag==RunMode);
            NumSolsNeed = 2*sum(vFixed(nO+1:nO+nD)==0);

            if NumFeasibleSols ~= NumSolsNeed
                warning('%d Solves, %d feasible solutions, but need %d. Terminating',NumSolvs,NumFeasibleSols,NumSolsNeed)
                return
            end

            %Prepare results for later use
            mResultsValUse = mResultsVal(vReturnFlag==RunMode,:);
            for i=1:n
                mResultsPlot(:,i) = ConvertPlot2DataUnits(mResultsValUse(:,i),mTransformToOrig(:,i),1);
            end
            mResCompact = [[min(vObjs(vReturnFlag==RunMode,:),[],1); max(vObjs(vReturnFlag==RunMode,:),[],1)] [min(mResultsValUse,[],1);max(mResultsValUse,[],1)]];
            
            if RunMode==2
                ColInds = ColInds(vFixed==0);
            end
        end  
        
        if lShowFiltered
            mResultsPlot = [mResultsCurrRec;mResultsPlot];
            lWidth = 1.5;
        else
            mResultsPlot = mResultsCurrRec;
        end
        
        %Show darkened lines for the current identified solutions (sub-set)
         if ~isempty(mResultsPlot) && ~any(any(isnan(mResultsPlot)))
             
            %get(hWindReturn,'OuterPosition')
            %get(ax2)
            SubSpacePCPDec = parallelcoords(mResultsPlot,'color',mColors(gDexes(1),1,3,:),'LineStyle',vLineStyleToUse,'linewidth',lWidth,'Standardize',UseStandarize);
            %get(ax2)
            
            if (nO>0)%%Overplot objectives in a different color
                SubSpacePCPObj = parallelcoords(mResultsPlot(:,1:nO+1),'color',mColors(gDexes(1),2,3,:),'LineStyle',vLineStyleToUse,'linewidth',lWidth,'Standardize',UseStandarize,'Visible',sObsVis);
            end
         else
             SubSpacePCPDec = parallelcoords([],'color',mColors(gDexes(1),1,2,:));
             if (nO>0) %%Overplot objectives in a different color
                SubSpacePCPObj = parallelcoords([],'color',mColors(gDexes(1),2,2,:));
             end
         end
         
         %Update the sliders
         NewLimitsPU = zeros(2,n);
         [vSliderValuesDU vSliderSteps vCurrValsPU vMins vMaxes] = ReadSliderValues(sSlider,vFixedVals,mTransformToOrig);
         
         NewLimitsPU = [vMins';vMaxes'];
         %[size(vSliderValuesDU);size(vSliderSteps);size(vCurrValsPU);size(NewLimits)]
         %vCurrValsPU = zeros(1,n);

         if (lHideSliderValue ~= 1) && ~isempty(mResCompact) && ~any(any(isnan(mResCompact)))
            for i=1:n
                for j=1:2
                    NewLimitsPU(j,i) = ConvertPlot2DataUnits(mResCompact(j,i),mTransformToOrig(:,i),1);
                end

                %vCurrValsPU(i) = ConvertPlot2DataUnits(vFixedVals(i),mTransformToOrig(:,i),1);
            end   
         end
         
         [sSliderOut,vShowSlider] = RenderSliders(NewLimitsPU(1,:),NewLimitsPU(2,:),vCurrValsPU',vSliderStep,sSliderUse);
                  
         %vPrint = [[1:n]' mResCompact(1,:)' vFixedVals mResCompact(2,:)']
    end

    function SetSliderValue(hind,event,sSlider,mConvertFactors,i)
        %called when a txtSliderValue value is entered
        %If eimSetBeyond is unchecked, checks entry is allowed (within slider range) and then sets the corresponding sSlider value
        %i is the index number of the txtSliderValue and sSlider
        %
        %Remember, txtValues are in original Data units and sSlider values are
        %in Plot units so we use mCovertFactors to convert

        %lBeyondExtents = get(cbAllowSets,'Value');
        lBeyondExtents = Enabled2Boolean(get(uimSetBeyond,'Checked'));
        
        sValueDU = str2double(get(hind,'String'));
        sValuePU = ConvertPlot2DataUnits(sValueDU,mConvertFactors(:,i),1);
        vFixedVals(i) = sValueDU;
        
        %Check within bounds
        sMaxValuePU = get(sSlider(i),'Max');
        sMinValuePU = get(sSlider(i),'Min');   
         
        if lBeyondExtents==0
            sDir = 'Nothing';
            
            if (sValuePU < sMinValuePU) 
                sValuePU = sMinValuePU;
                sDir = 'minimum';
            end
            if (sValuePU>sMaxValuePU)
                sValuePU = sMaxValuePU;
                sDir = 'maximum';
            end

            sValuesPU = [sValuePU sMinValuePU sMaxValuePU];
            sValuesDU = ConvertPlot2DataUnits(sValuesPU,mConvertFactors(:,i));

            vFixedVals(i) = sValuesDU(1);
            
            if sDir ~= 'Nothing'
                %outside of range
                %h = errordlg(sprintf('Input (%.2f) is outside allowable range: %.2f to %.2f. Setting to %s value (%.2f). Please change if needed.',sValueDU, sMinValue,sMaxValue,sDir,sSetValue),'Input Error');
                h = errordlg(sprintf('Input (%.2f) is outside allowable range: %.2f to %.2f. Setting to %s value (%.2f). Please change if needed.',sValueDU, sValuesDU(2),sValuesDU(3),sDir,sValuesDU(1)),'Input Error');

                uiwait(h);
                %return to extent value
                set(sSlider(i),'Value',sValuePU);
                set(hind,'String',sprintf('%.2f',sValuesDU(1)));

                uicontrol(hind);
                return
            end
            
            set(sSlider(i),'Value',sValuePU);
            UpdateHighlightedTraces
        else
            if sValuePU < sMinValuePU
                sMinValuePU = sValuePU;
            else
                sMaxValuePU = sValuePU;
            end
            [sSlider(i),vShowSliders(i)] = RenderSliders(sMinValuePU,sMaxValuePU,sValuePU,vSliderStep(i),sSlider(i),i-1,length(sSlider));
        end
    end

    function SetTxtValues(hind,event,txtSliderValue,mConvertFactors,i)
        %called when a sSliderValue value is changed
        %and then sets the corresponding txtSliderValue value
        %i is the index number of the txtSliderValue and sSlider

        %Remember, txtValues are in original Data units and sSlider values are
        %in Plot units so we use mCovertFactors to convert


        %txtSliderValue;

        sValuePU = get(hind,'Value');
        sValueDU = ConvertPlot2DataUnits(sValuePU,mConvertFactors(:,i));

        set(txtSliderValue(i),'String',sprintf('%.2f',sValueDU));
        
        lHideSliderValue = Enabled2Boolean(get(uimHideSliders,'Checked'));
        
        if lHideSliderValue==0
            %also set the fixed checkbox for the axis
            set(cbChecks(i),'value',1);
        end
        
        vFixedVals(i) = sValueDU;       
        UpdateHighlightedTraces        
    end   

    function AxisChecked(hind,event,cbCheck,i)
        %Callback when an axis is checked/unchecked. Fires the update
        UpdateHighlightedTraces
    end

    function callHideSliders(hind,event,sSlider)
        %hides all sliders when checked, makes them visible when unchecked
        sValue = ToggleMenuItem(hind,0);
        
        if sValue==1
            sVisible='off';
        else
            sVisible='on';
        end

        [m n] = size(sSlider);

        for i=1:n
            if sSlider(i) > 0
                set(sSlider(i),'Visible',sVisible);
            end
        end
    end

    function ToggleCheckBoxVisibility(hind,event,cbChecks,hTexts,hTabContainer)
        %when hind is checked, hides all the checkboxes below the axes and moves the
        %axis text labels up to the yPosChecked value.

        %When unchecked, makes the checkboxes visible and returns the text
        %labels to the yPosUnchecked value

        %Similarly shifts the GroupTable
        
        % INPUTS
        %  hind = handle to the menu item that controls this function
        %  cbChecks = vector of handles to the checkbox controls
        %  hTexts = vector of handles to the text labels thta sit under the
        %       checkboxes
        %  hTabContainer = handle of the container which holds the table
        %       for grouping axes labels
        
        %If sliders are on, can not hide checkboxes      
        if (Enabled2Boolean(get(hind,'Checked'))==0) &&  (Enabled2Boolean(get(uimHideSliders,'Checked'))==0)
            %Turn on the checks
            sResponse = questdlg('Sliders are visible. Can not hide checkboxes. Hide both?','How to proceed','Yes','Cancel','Cancel');
            
            if strcmpi(sResponse,'Yes')
                %turn off the slideers
                ToggleCurrMenuItem(uimHideSliders,0);
            else %don't do anything
                return
            end
        end

        sValue = ToggleMenuItem(hind,0);%get(hind,'Value');
        set(cbChecks(1),'Units','Characters');
        vChecksPos = get(cbChecks(1),'Position');
        set(cbChecks(1),'Units','Normalized');
        
        if sValue==1
            sVisible='off';
            yChange = vChecksPos(4);
        else
            sVisible='on';
            yChange = -vChecksPos(4);
        end

        [m n] = size(cbChecks);
        
        %Hide the checks, reposition the texts
        for i=1:n
            set(cbChecks(i),'Visible',sVisible);
            cUnits = get(hTexts(i),'Units');
            set(hTexts(i),'Units','characters');
            cPos = get(hTexts(i),'Position');
            cPos(2) = cPos(2) + yChange;
            set(hTexts(i),'Position',cPos);
            set(hTexts(i),'units',cUnits);
        end

        %Reposition the TabContainer
        %hParent = get(hTabContainer,'Parent');
        tUnits = get(hTabContainer,'Units');
        set(hTabContainer,'Units','Characters');
        vPos = get(hTabContainer,'Position');
        vPos(2) = vPos(2) + yChange;
        set(hTabContainer,'Position',vPos);
        set(hTabContainer,'Units',tUnits);
    end

    function ShowAllControls(hind,event,cFrame,hAxes,cChecks,sSlides,hTabContainer,xWidthAddWOControls)
        %Toggles the visibility of all the controls
        %if off, then control frame is off
        %if on, then main on and slider frame returned to prior value
        %also adjusts the width of the main plot and the positions of the
        %checkboxes and slider controls (on the graph)

        vPos = get(hAxes(1),'Position');
        xLeft = vPos(1);
        yBottom = vPos(2);
        xWidth = vPos(3);
        yHeight = vPos(4);

        nAxes = max(size(hAxes));
        nControls = max(size(cChecks));
        nSlides = max(size(sSlides));
        %nLines = max(size(hLines));

        xWidthOld = xWidth;

        %Turn the uimenu (hind) to the opposite of what it currently is
        strVis = get(hind,'Checked');
        cCheckedSet = Boolean2Enabled(ToggleMenuItem(hind,event));

        if strcmpi(cCheckedSet,'on') %cVal==0 %make all the frames invisible                 
            xWidth = xWidth + xWidthAddWOControls; 
        else
            xWidth = xWidth - xWidthAddWOControls;
        end

        set(cFrame,'Visible',strVis);

        %set visibility on all the children panels and their controls
        SetVisibilityAll(cFrame,strVis);

        %position the plot areas
        for i=1:nAxes
            set(hAxes(i),'Position',[xLeft yBottom xWidth yHeight]);
        end

        %position the checkbox and slider controls (on the graph)
        vCheckPos = get(cChecks(1),'Position');

        for i=1:nControls
            vCheckPos(1)=xLeft+(YAxisMargins(1)+i-1)*xWidth/(nControls-1+sum(YAxisMargins))-0.005;

            set(cChecks(i),'Position',vCheckPos);
            if sSlides(i)~=0
                vSliderPos = get(sSlides(i),'Position');
                vSliderPos(1)=xLeft+(YAxisMargins(1)+i-1)*xWidth/(nSlides-1+sum(YAxisMargins))-0.005;
                set(sSlides(i),'Position',vSliderPos);
            end
        end

        PositionGroupLabelsTable(hTabContainer, nO, hTexts, hAxes(1));
        UpdateLegend;
    end

    function HideCheckboxesOld(hind,event,cbChecks,hTexts,hTabContainer,yPosUnchecked,yPosChecked,yPosGroupOffset, yPosGroupOffsetUnits)
        %when checked, hides all the checkboxes below the axes and moves the
        %axis text labels up to the yPosChecked value.

        %When unchecked, makes the checkboxes visible and returns the text
        %labels to the yPosUnchecked value

        % yPosGroupOffset -- the measure in yPostGroupOffsetUnits by which to
        % vertically shift the table of group text labels

        sValue = ToggleMenuItem(hind,0);%get(hind,'Value');

        if sValue==1
            sVisible='off';
            yPos = yPosChecked;
            %shift the table up by yPosGroupOffset
            yPosGroupChange = yPosGroupOffset;
        else
            sVisible='on';
            yPos = yPosUnchecked;
            yPosGroupChange = -yPosGroupOffset;
        end

        [m n] = size(cbChecks);
 
        for i=1:n
            set(cbChecks(i),'Visible',sVisible);
            cUnits = get(hTexts(i),'Units');
            set(hTexts(i),'Units','data');
            cPos = get(hTexts(i),'Position');
            set(hTexts(i),'Position',[cPos(1) yPos]);
            set(hTexts(i),'units',cUnits);
        end


        %hParent = get(hTabContainer,'Parent');
        tUnits = get(hTabContainer,'Units');
        vPos = get(hTabContainer,'Position');

        if ~strcmp(tUnits,yPosGroupOffsetUnits)
            %need to convert to correct units
            set(hTabContainer,'Units',yPosGroupOffsetUnits);
            vAlt = get(hTabContainer,'Position');
            yPosGroupOffset = yPostGroupOffset*vPos(2)/vAlt(2);
            %return the units back
            set(hTabContainer,'Units',pUnits);
        end

        vPos(2) = vPos(2)+yPosGroupChange;
        set(hTabContainer,'Position',vPos);

        %for i=1:nG
        %    cPos = get(hTextGroup(i),'Position');
        %    set(hTextGroup(i),'Position',[cPos(1) cPos(2)+yPosGroupChange]);
        %end

        %yPosGroupChange;
        %vChanges = ConvertAxisToFig([0 yPosGroupChange],get(hTextGroup(1),'Parent'));

        %for i=1:nL
        %    if hLines(i) > 0
        %        cPos = get(hLines(i),'Position');
        %        set(hLines(i),'Position',[cPos(1) vChanges(2) cPos(3) cPos(4)]);
        %    end
        %end
    end

    function CheckAllBoxes(hind,event,cbTargetChecks)
        %called when user checks the cbCheckAll box
        %sets all the axes checkboxes cbChecks to the same value

        CheckValue = ToggleMenuItem(hind);
            
        [m n] = size(cbTargetChecks);
        for i=1:n
            set(cbTargetChecks(i),'Value',CheckValue);
        end
        UpdateHighlightedTraces
    end

    function SetAxesInFilterAltsBox(hind,event,txtSlides,lblSlides,nO,startHeightPx,hPerRow)
       % Handler for when a new values is entered for the 1st Decision Axis to List on the
       % Filter Alternatives menus. Shows the specified axis in first and
       % all below until run to bottom of the pane.
       
       try
          lTry = max([floor(str2num(get(hind,'string'))) 1]);
       catch
          msgbox('Error with input %s. Must be numeric',get(hind,'string'));
          return
       end
       
       nFull = length(txtSlides);
       
       if lTry>(nFull-nO)
           lStart = nFull;
       else
           lStart = nO+lTry;
       end
       
       for i=nO+1:nFull
           sVis = 'on';
           vPosVert = startHeightPx - hPerRow*(i-lStart);
           if (i<lStart) || (vPosVert < 0)
              sVis = 'off';
           end
           tUnits = get(txtSlides(i),'units');
           set(txtSlides(i),'units','pixel');
           vTxtPos = get(txtSlides(i),'Position');
           vTxtPos(2) = vPosVert;
           lUnits = get(lblSlides(i),'units');
           set(lblSlides(i),'units','pixel');           
           vLblPos = get(lblSlides(i),'Position');
           vLblPos(2) = vPosVert;
           
           set(txtSlides(i),'Visible',sVis,'Position',vTxtPos);
           set(lblSlides(i),'Visible',sVis,'Position',vLblPos);
           
           set(txtSlides(i),'units',tUnits);
           set(lblSlides(i),'units',lUnits);
       end       
    end

    function UpdateLegendMenu(hind,event,i,Total)
        % Handler for when a legend menu option is selected.
        % Makes sure only the current menu option i is selected
        % Then calls the update legend handler
        
        set(uimLegendOptions([1:Total] ~= i),'Checked','off');
        set(uimLegendOptions(i),'Checked','on');
        UpdateLegend;
    end

    function UpdateLegend
        % Shows / hides the legend depending on the settings for
        %   - ShowLegend
        %   - ShowControls
        %   - ShowPareto
        % Also sets the entries in the legend from the groups checked
        % (showing)
        
        [vParams,hControls] = ReadControls('UpdateLegend',hWindReturn);
        [cShowLegend,cShowControls,cShowInset, cmGroupData] = aGetField(vParams,{'ShowLegend','ShowControls','ShowInsetPlot', 'mGroupData'});
        
        lnU = size(cmGroupData,1);
        %Construct a handle to the lengend objects
        hPCTraces = zeros(lnU,1);
        for i=1:lnU
            hPCTraces(i) = hPCGroup{i,1}(1);
        end
        checkedGroups = cell2mat(cmGroupData(:,2));
        
        % Set visibility
        visSet = 'show';
        if (cShowInset==1) || strcmpi(cShowLegend,cLegendOptions{1}) || ...
                (strcmpi(cShowLegend,cLegendOptions{2}) && (cShowControls==1)) || ...
                (get(cbRampColor,'value') ==1)
            %Hide if inset plot is visible OR Show Legend is set to hide OR
            %    Show Legend is set to only when control panel hidden and
            %    control panel is hidden OR Color Ramp box is checked
            visSet = 'hide';
        else
            hLegend = legend(hPCTraces(checkedGroups==1),cmGroupData(checkedGroups==1,1),'FontSize',FontSize-4);
        end
        if hLegend > 0
            legend(hLegend,visSet);
        end
    end
            
    function ToggleSliders(hind,event)
        % Handler for when the slider menu item is selected
        % If going to show sliders, we also need to show the axis checkboxes
        if (Enabled2Boolean(get(hind,'Checked'))==1) &&  (Enabled2Boolean(get(uimHideChecks,'Checked'))==1)
            %Turn on the checks
            ToggleCheckBoxVisibility(uimHideChecks,0,cbChecks,hTexts,hTabContainer);
        end
        
        ToggleCurrMenuItem(hind,event);
    end

    function ToggleCurrMenuItem(hind,event)
        % toggles whether to show the trace for the current selected record
        % and highlighted records
        cVal = ToggleMenuItem(hind); %cVal = get(hind,'value');
        
        UpdateHighlightedTraces;
    end

   function ToggleShowObjs(hind,event)
        % toggles whether to show the objective function axes in a
        % different color
        cVal = ToggleMenuItem(hind); %cVal = get(hind,'value');
        ShowGroups(0,0,hWindReturn,2);
   end

    function GroupDecisionLabels(hind,event,hTexts,hTabContainer, nO, hAxis)
        % toggles the turning on/off the hTabContainer which contains the
        % axis labels in a grouped html table
        %
        % hTexts are the handles to the individual axis labels that appear
        %   below the axis
        % hTabContainer is the handle to the java element showing the table of
        %   grouped decision variables
        % nO is an offset for the first nO elements of hTexts to position the left side the table container 
        % hAxis is handle to the axis to which the Text Labels are
        %   associated

        % also recalculates the position for the hTabContainer table just below
        % the hTexts

        cVal = ToggleMenuItem(hind); %cVal = get(hind,'value');
        gText = Boolean2Enabled(cVal);
        
        %Recalculate and set the tables position
        PositionGroupLabelsTable(hTabContainer,nO,hTexts,hAxis);
       
        set(hTabContainer,'Visible',gText);
    end

    function [vTablePosition] = PositionGroupLabelsTable(hTabContainer, nO, hTexts, hAxis)
        % Position the hTabContainer based on the current location of the
        % axis and text labels. This position is:
        %   Left: nO+1 axis
        %   Bottom: Below the bottom of the hTexts
        %
        % hTexts are the handles to the individual axis labels that appear
        %   below the axis
        % hTabContainer is the handle to the java element showing the table of
        %   grouped decision variables
        % nO is an offset for the first nO elements of hTexts to position the left side the table container 
        % hAxis is handle to the axis to which the Text Labels are
        %   associated
        % OUTPUTS
        %   vTablePosition = the newly calculated position vector for the hTabContainer

        % Get the extents of the texts
        mTextExtent = cell2mat(get(hTexts(nO+1:end),'Extent'));
        mMaxExtent = min(mTextExtent(:,2)); %Read the minimum bottom extent (lowest)
        %Remember, text extents are relative to the axis, the Table is
        %relative to the figure!
        
        %Now work from the axis position
        aPosition = get(hAxis,'Position');
        xLeftTab = aPosition(1)+(YAxisMargins(1)+(nO-1)+0.5-.0025)/(n-1+sum(YAxisMargins))*aPosition(3);
        xWidthTab = aPosition(3)*(nD)/(n-1+sum(YAxisMargins));
        vTabPosCur = get(hTabContainer,'Position');
        TabFontSizes = FontSize-4-12+1+[1:lNumFields];
        lTabHeight = sum(TabFontSizes/TabFontSizes(1));
        yHeightTab = lTabHeight*0.06; %0.07 if a row wraps in the table
      
        yBottomTab = max([aPosition(2)-aPosition(4)*(-mMaxExtent) - yHeightTab] - 0.01);      %0.025     
        %position the table below the largest axis text label
        vTablePosition = [xLeftTab yBottomTab xWidthTab yHeightTab];
        set(hTabContainer,'Position',vTablePosition);
    end

    function [vSlidersOut,vShowSlider] = RenderSliders(vMins,vMaxes,vFixedValsPU,vSliderStep,vSliders,Indexes,TotalInds) 
        % Renders the height of the sliders according to the specified
        % range in vMins (minimums) and vMaxes(maximumns). Hides sliders that have the same min and max
        %
        % INPUTS
        % vFixedValsPU = is the current setting of the slider in Plot Units
        % vSliderStep = is the step for the slider to move when clicking an
        %   end (in plot units)
        % Optional vSliders is an array of object ids to the sliders. If
        %   passed, the function will update. If omited, the function will
        %   create new arrays
        % Optional Indexes are the original indexes of the sliders (if a
        %   subset are passed) to position correctly
        % Optional TotInds is the total number original sliders (if a
        % subset are passed)
        % 
        
        % OUTPUTS
        % vSlidersOut = vector of object IDs for the sliders (same as
        % vSliders, if  vSliders is passed
        % vShowSlider = vector of binary (1=show slider, 0=hide slider).
        %       Hide if the min and max for the decision variable are the same
        
        % Find the size of the vectors
        n = size(vMins,2);
        
        % Echo the inputs
        %[size(vMins); size(vMaxes); size(vFixedValsPU); size(vSliderStep)]
        %if nargin==5
        %    size(vSliders)
        %end
        %[vMins; vMaxes; vFixedValsPU;vSliderStep]
        
        if (n~=size(vMaxes,2)) || (n~=size(vFixedValsPU,2)) || (n~=size(vSliderStep,2)) || ((nargin==5) && (n~=size(vSliders,2)))
            error('RenderSliders: Paremeters are incompatible sizes, vMins: %d %d, vMaxes: %d %d, vVixedValsPU: %d %d, vSliderStep: %d %d',size(vMins),size(vMaxes),size(vFixedValsPU),size(vSliderStep))
            return
        end
        
        %Initialize variables
        vShowSlider = vMins < vMaxes;
        
        lHideSliderValue = Enabled2Boolean(get(uimHideSliders,'Checked')); %get(cbHideSliders,'value');
        
        if nargin>=5
            vSlidersOut = vSliders;
        else
            vSlidersOut = zeros(1,n);
        end
        
        if nargin < 6
           Indexes = zeros(n,1); 
           TotalInds = n;
        end
        
        %Loop through the sliders
        for i=1:n
            vSliderStepValue = [vSliderStep(i) vSliderStep(i)];

            sVis = 'on';
            Delta = 0;
            
            if nargin < 5
                vSlidersOut(i) = uicontrol('Style','slider');
            end
            
            if (lHideSliderValue==1) || (nargin<5) 
                %The checkbox value overrules everything else. Also, don't
                %show if initially creating
                sVis = 'off';
            end
            
            %Check the fixed value is within the limits. If not, set to the
            %nearest limit
            SliderSetValue = vFixedValsPU(i);
             
            if vMaxes(i) < vMins(i)
                vMaxes(i) = vMins(i);
            end
            
            if (vMins(i) > SliderSetValue) 
                  SliderSetValue = vMins(i);
            end
            
            if (vMaxes(i) < SliderSetValue) 
                  SliderSetValue = vMaxes(i);
            end
            
            if (vShowSlider(i)==0) || (vMaxes(i)==vMins(i))
                sVis = 'off';
            end
            
            if vMaxes(i) == vMins(i)
                Delta = 0.01;
            end
            
            %[i Delta xLeft+(i-1)*xWidth/(n-0.5)-0.005 yBottom+yHeight*(vMins(i)-ymin)/(ymax-ymin) 0.01 yHeight*(vMaxes(i)+Delta-vMins(i))/(ymax-ymin)]
            
            set(vSlidersOut(i),'Min',vMins(i),'Max',vMaxes(i)+Delta,'Value',SliderSetValue,'SliderStep',vSliderStepValue,'Visible',sVis,...
                    'Units','normalized','Position', [xLeft+(YAxisMargins(1)+i+Indexes(i)-1)*xWidth/(TotalInds-1+sum(YAxisMargins))-0.005 yBottom+yHeight*(vMins(i)-ymin)/(ymax-ymin) 0.01 yHeight*(vMaxes(i)+Delta-vMins(i))/(ymax-ymin)]);
        end    

    end

    function SetSliderValsToRecord(hind,event,RecordInd)
         %sets all the sliders and associated text boxes to the vector of
         %values in vSetValsPU
         %
         % vSetValsPU = vector of values in Plot Units to set.

         %Remember, sliders are in Plot Units and txtTexts are in Data units.
         %Use mConvertFactors to convert from Plot to Data units.

         n = length(sSlider);    

         %Update the record control
         set(txtCurrRec,'string',num2str(RecordInd));

         %fprintf('Before\n');
         %[vFixedValCurr vSliderSteps vFixedValPU] = ReadSliderValues(sSlider,zeros(1,n),mTransformToOrig);

         for i=1:n
            if vShowSliders(i)>0       
                vValue = mPlot(RecordInd,i);
                
                set(sSlider(i),'Value',vValue);
                vValueDU = ConvertPlot2DataUnits(vValue,mTransformToOrig(:,i));
                set(txtSliderValue(i),'String',sprintf('%.2f',vValueDU));            
                vFixedVals(i) = vValueDU;
            end
         end
         UpdateHighlightedTraces
    end

    function TriggerUpdateHighlightTraces(hind,event)
        UpdateHighlightedTraces
    end

    function RenameGroup(hind,event,i)
        % Takes the newly renamed i'th element in txtGroupNames text box
        % and changes it in vGroup and vUnique
        NewName = get(txtGroupNames(i),'string');
        
        %Find the rows in vGroup that have the old name        
        if iscell(vGroup)
            GrpRows = strcmpi(vGroup,vUnique(i));
        else
            GrpRows = vGroup==vUnique(i);
        end
        lNumRows = sum(GrpRows);
        
        %Replace those rows with the new name        
        vGroup(GrpRows) = repmat(NewName,lNumRows,1);
        %Update vGroup in the Figure's application data stream
        vararg_curr = aSetField(vararg_curr,'vGroup',vGroup);
        setappdata(hWindReturn,'varargs',vararg_curr);
        %Set the local variable for Group Names to the new name
        vUnique(i) = NewName;        
    end

   function [hPCPs,cbColors] = CreateEmptyColorRampPCPs(txtNumColors,lNumColors)
        % Creates a cell array of pointers to n empty parallel
        % coordinate plots where n is the number of color classes specified
        % in txtNumColors. Also checkboxes on the Display tab
        % to toggle each color class. lNumColors is optional parameter with
        % number of classes specified as a double
        
        % OUTPUTS
        %   hPCPs = cell array of pointers
        %   cbColors = array of handles to the color check boxes
       
       %Read in the number of color classes
       
       if nargin == 2
           lNumClasses = lNumColors;
       else
           lNumClasses = round(str2num(get(txtNumColors,'string')));

           if isempty(lNumClasses)
               warning('Invalid input for # Color Classes (%s). Defaulting to 10',get(txtNumColors,'string'));
               lNumClasses=10;
           end
       end
       
       lNumClasses = round(lNumClasses);
       
       hPCPs = cell(lNumClasses,1);
       cbColors = zeros(lNumClasses,1);
       hPCPs(:) = {0};
       
       
       %for i=1:lNumClasses
       %    hPCPs{i} = parallelcoords([]);
      %     cbColors(i) = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', ['Class ',num2str(i)],'Position',[25 lTopColor-25-(i)*20 230 20],'fontsize', fontsizecntls-2,'Callback',{@ToggleColor,i},'Value',1,'visible','off');
      % end      
    end

    function [hPCPout,cbColorOut] = DeleteColorRampPCPs(hPCPIn,cbColors,txtNumClasses)
        %deletes the parallel coordinate objects in hPCPIn
        
        if (length(hPCPIn) > 0) && (length(cbColors)>0)
            for i=1:length(hPCPIn)
                if hPCPIn{i} > 0
                    delete(hPCPIn{i});
                end
                if cbColors(i) > 0 
                    delete(cbColors(i));
                end
            end
        end
        
        %Reseed
        [hPCPout,cbColorOut] = CreateEmptyColorRampPCPs(txtNumClasses);
    end

    function ToggleColor(hind,event,i)
        %Toggle on/off the PCP color class associated with the i-th checkbox

        lValue = get(hColorChecks(i),'Value');
        
        if lValue == 1
            %Turn it on, make it visible
            sVis = 'on';
        else
            sVis = 'off';
        end
        
        set(ColorRampPCPs{i},'visible',sVis);       
    end

    function ShowTab(hObj,event,i,hTabs,hTabButtons,vBackColor)
       %selects and shows the ith member of hTabs and hTabButtons
       %turns the color of the current button to white and the other buttons to
       %the specified background color
       %
       % i = index of the tab selected to show
       % hTabs = vector of handles to the tabs
       % hTabButtons = vector of handles to the tab buttons
       % vBackColor = vector specifying color to turn background of
       % non-selected buttons

       nTabs = max(size(hTabs));

       for k=1:nTabs
           if i==k
               sVis = 'on';
               sColor = 'white';
           else
               sVis = 'off';
               sColor = vBackColor;
           end

           %k;
           %sColor;

           set(hTabs(k),'visible',sVis);
           set(hTabButtons(k),'BackgroundColor',sColor);
       end
       if i==2 && (Enabled2Boolean(get(uimHideChecks,'checked'))==1)
           %Additionally turn on the axis check boxes
           ToggleCheckBoxVisibility(uimHideChecks,0,cbChecks,hTexts,hTabContainer);
       end
    end


    function RampColor(hind,event)
       % If the check box is checked, ramps the colors of traces in the specified direction based on
       % values on the first checked axis by turning on a new layer
       % Otherwise, turns off the ramping layer
       
       % Direction => Ascend: lightest to darkest (smallest value to
       %                largest value; darkest layer on top/last)
       %              Descend: reverse -- darkest to lightest (largest
       %              value to smallest value; darkest layer on top/last)
       
       % Read the color ramp tab control settings
       lChecked = get(cbRampColor,'value');
       lDirection = get(cbDirection,'value');
       lHighlightGroup = get(cbHighlightGroup,'Value') - 1;
       
       % Delete the existing color ramp layer
       [ColorRampPCPs, hColorChecks] = DeleteColorRampPCPs(ColorRampPCPs,hColorChecks,txtNumClasses);
       
       hold on
       
       if lChecked == 0
           %turn off the color ramp layer and turn back on the original
           %layers
           ShowGroups(0,0,hWindReturn,2);
       else
           %turn off the prior layers and turn on the color ramp layer for the first checked axis
           
           %turn off the regular plot layers
           for l = 1:size(hPCGroup,1)
               if 1 %l ~= lHighlightGroup % Hide all groups except the highlight group
                   for j=1:2
                       for k=1:max(size(hPCGroup{l,j}))
                           set(hPCGroup{l,j}(k),'Visible','off');
                       end
                   end
               end
           end
           
           %Find the Axis and number of color classes
           %[sFirst,lFirst] = GetFirst(cbChecks, {vObjLabels{:} vXLabels{:}}, 1);
           [sFirst,lFirst] = GetFirst(cbChecks, get(hTexts,'string'), 1);

           %Read in the number of color classes
           lNumClasses = str2num(get(txtNumClasses,'string'));

           if isempty(lNumClasses)
               warning('Invalid input for # Color Classes (%s). Defaulting to 10',get(txtNumClasses,'string'));
               lNumClasses=10;
           end
           
           %Figure which rows in the matrix to use from the groups seleted
           [vParams,hControls] = ReadControls('ShowGroups',hWindReturn);
           [mGroupInfo,vGroupCurr] = aGetField(vParams,{'mGroupData' 'vGroup'});
           GroupsToUse = mGroupInfo((cell2mat(mGroupInfo(:,2))==1),1);         
           vRowToUse = zeros(length(vGroupCurr),1);
    
            for i=1:size(GroupsToUse,1)
                vRowToUse = vRowToUse+strcmpi(vGroupCurr,GroupsToUse{i});
            end
                     
           ranges = [min(mPlot(vRowToUse==1,lFirst)) max(mPlot(vRowToUse==1,lFirst))];
               
           if ranges(1) == ranges(2)
               warning('RampColor: Min and Max on %s (%d-th axis) are the same. Can not ramp color. Select a different axis.',sFirst{:},lFirst)
               %turn off the color ramp layer      
               return
           end

           lNumClasses = round(lNumClasses);   

           %Build the color ramp with the specified number of classes          
           %For blues
           %BaseRamp = (OSUColorRamps('LightToDarkBlue7Step')); %reshape(mColors(1,1,:,1:end),3,3)
           %BaseRamp = BaseRamp(2:7,:);
           %For Green
           BaseRamp = (OSUColorRamps('GreenToMagenta16Step'));
           BaseRamp = flipud(BaseRamp(1:7,:)) ;
           
           mRamp = ExpandColorRamp(lNumClasses,BaseRamp,0,0);
           
           %Build the position matrix
           
           vOffset = mod([1:lNumClasses]-1,lNumRowsForColorChecks-1)+1;
           hOffset = floor(([1:lNumClasses]-1)/(lNumRowsForColorChecks-1));
                      
           if lDirection == 2
               %descend
               cRange = [lNumClasses 1];
               iAdd = 1;
           else
               %ascend
               cRange = [1 lNumClasses];
               vOffset = fliplr(vOffset);
               hOffset = fliplr(hOffset);
               iAdd = 0;
           end

           %Reasign the data to the classes. Use linear interpolation along the
           %axis
           vGroupNew = round(interp1(ranges,cRange,mPlot(:,lFirst)));
           DataRanges = [min(mResults(vRowToUse==1,lFirst)) max(mResults(vRowToUse==1,lFirst))];
           
           ClassSignifMag = floor(log((DataRanges(2) - DataRanges(1))/lNumClasses)/log(10)); %Magnitude of the difference between class values
           
           vGroupRanges = interp1(cRange,DataRanges,[1:(lNumClasses-1)/lNumClasses:lNumClasses]);
           vGroupPercents = (vGroupRanges)/vGroupRanges(lNumClasses+1)*100;
           
           %Cycle through the classes, create one plot and checkbox per color/group
           lWidth = 250;
           for i=1:lNumClasses
                ColorRampPCPs{i} = parallelcoords(plot2,mPlot(((vGroupNew==i) + (vRowToUse==1)==2),:),'color',mRamp(i,:));
                %Reset the parent because parallelcoords doesn't seem to do
                %so
                set(ColorRampPCPs{i}(:),'Parent',plot2);
                vPosition = [10+(lWidth+5)*hOffset(i) lTopColor-65-(vOffset(i))*20 lWidth 20];
 
                %As actual value
                sStartVal = ThousandSep(10^ClassSignifMag*round(vGroupRanges(i+iAdd)/10^ClassSignifMag));
                sEndVal = ThousandSep(10^ClassSignifMag*round(vGroupRanges(i+1-iAdd)/10^ClassSignifMag));
                
                %As percent
                sStartPer = (vGroupPercents(i+iAdd));
                sEndPer = (vGroupPercents(i+1-iAdd));
                
                hColorChecks(i) = uicontrol('Parent',hTabs(2),'Style', 'checkbox', 'Position',vPosition, 'ForegroundColor',[0 0 0],...
                    'fontsize', fontsizecntls,'Callback',{@ToggleColor,i},'Value',1,'visible','on', ...
                    'string',sprintf('%s to %s (%.f%% to %.f%%)',sStartVal{:},sEndVal{:},sStartPer,sEndPer),'BackgroundColor',mRamp(i,:));
           end
           
           if lHighlightGroup > 0
               drawnow; pause(0.1); %let the gui catch up
               uistack(cell2mat(hPCGroup(lHighlightGroup,1)),'top');
           end
           UpdateLegend;
       end       
    end

    function [lIndToUse,lLineWidth,mColorsDecs,mColorsObjs,strVis,sMarker] = GetGroupRenderings(lGroupNumber,cMarkers)
        %Returns the index, line width, decision and objective function
        %color, and visibility needed to render the lGroupNumber group
        %Also the marker if the optional cMarkers is passed
        
        lIndToUse = 1+mod(lGroupNumber-1,size(mColors,1)); %vUnique(i)

        if lOptGroupCurr == lGroupNumber                
            lLineWidth=2;
            vGroupInds = [1:nU]';
            lLineWidth = max([vGroupThick(lOptGroupCurr) max(vGroupThick(vGroupInds~=lGroupNumber))+1]); %one more than the larger thickness that is not the highlightgroup
            mColorsDecs = mHighlightColor;
            mColorsObjs = mHighlightColor;
            sMarker = '^';
        else
            lLineWidth=vGroupThick(i);
            mColorsDecs = mColorsDecsGroup(lIndToUse,:);
            mColorsObjs = mColorsObjsGroup(lIndToUse,:);
            if ~isempty(cMarkers)
                sMarker = cMarkers{i};
            end
        end

        strVis = Boolean2Enabled(vShowGroup(lGroupNumber));       
    end

    function ShowGroups(hObj,event,hWindCurr,HighlightThicknesses)
        % Shows the traces for the groups that are checked; hides traces for
        %       groups that are not checked

        % hWindCurr = Handle to the current Figure
        % HighlightThickness = line thickness for group to highlight

        [vParams,hControls] = ReadControls('ShowGroups',hWindCurr);
        [mGroupInfo,lShowObjsDiffColor] = aGetField(vParams,{'mGroupData' 'ShowObjsDiffColor'});
        [cbGroupsOneColor] = aGetField(hControls,{'AllGroupsSameColor'});
        %read in the check box values        
        cAllSameColor = get(cbGroupsOneColor,'Value');
        vGroupChecked = cell2mat(mGroupInfo(:,2));
        vGroupThicks = cell2mat(mGroupInfo(:,3));
        
        nU = length(vGroupChecked);
        
        for i=1:nU
            if vGroupChecked(i)==1
                strFullVis = 'on';
            else
                strFullVis = 'off';
            end
            
            strObjVis = Boolean2Enabled(vGroupChecked(i) && lShowObjsDiffColor);

            if cAllSameColor==1
                vGroupToUse = 1;
            else
                vGroupToUse = i;
            end   

            mColorToUse = squeeze(mColors(vGroupToUse,1:2,1,:));
            lLineThick = vGroupThicks(i);

            if (i==lOptGroupCurr)
                for j=1:2
                    mColorToUse(j,:) = mHighlightColor;
                end
                %lLineThick = HighlightThicknesses(1);      
                vGroupInds = [1:nU]';
                lLineThick =  max([vGroupThicks(lOptGroupCurr) max(vGroupThicks(vGroupInds~=i))+1]); %one more than the larger thickness that is not the highlightgroup

            end

            %make the group visible
            for j=1:2
               for k=1:max(size(hPCGroup{i,j}))
                   if j==1
                       strVis = strFullVis;
                   else
                       strVis = strObjVis;
                   end
                   
                   set(hPCGroup{i,j}(k),'Visible',strVis,'Color',mColorToUse(j,:),'linewidth',lLineThick);
               end
            end
        end
        UpdateSetToControls(0,0,hWindCurr,0,0); 
        UpdateInset;
        UpdateLegend;
    end

    function UpdateSetToControls(hind,but,hWindCurr,direction,magnitude)
        % Updates the SetTo controls (record number and number of records)
        % based on the currently selected groups and the direction and
        % magnitude to advance records
        %
        % direction: -1 (left), +1 (right), 0 (jump right to specified record)
        % magnitude: 1 (advance one) or 2 (advance to far end)

        [vararg_curr,hControls] = ReadControls('UpdateSetToControls',hWindCurr);
        [vGroup,iOpt,mGroupData,vGroup,lCurrRecord] =  aGetField(vararg_curr,{'vGroup' 'GroupToHighlight' 'mGroupData' 'vGroup' 'CurrentRecord'});
        [GroupToHighBut] = aGetField(hControls,{'GroupToHightlightBut'});

        % Filter mGroupData by selected groups
        mCheckedGroups = mGroupData(cell2mat(mGroupData(:,2))==1,1);

        %Enable the GroupToHighlightButton if it's in the filtered set
        if ismember(iOpt,mCheckedGroups)
            set(GroupToHighBut,'enable','on');
        else
            set(GroupToHighBut,'enable','off');
        end

        if iscell(vGroup)
            vRecordsToUse = ismember(vGroup,mCheckedGroups);
        else
            vRecordsToUse = ismember(vGroup,cell2mat(mCheckedGroups));
        end

        m=length(vGroup); 
        inds = [1:m]';
        CurrValidRecords = inds(vRecordsToUse==1); %Indexes of records that are in checked groups
        
        %If current record belongs to a group that is no
        %longer shown, advance to nearest record of a group shown
        if vRecordsToUse(lCurrRecord) == 0
            lDistFromCurr = [-lCurrRecord:m-lCurrRecord]';
            [lNear,lNearInd] = min(abs(lDistFromCurr(vRecordsToUse==1)));
            lCurrRecord = CurrValidRecords(lNearInd);
        elseif ~isempty(lCurrRecord)
            lNearInd = find(CurrValidRecords==lCurrRecord);
        end

        if direction ~= 0 
            %Advance to the specified records
            if (magnitude == 1) && (((direction == -1) && (lNearInd > 1)) || ((direction == 1) && (lNearInd < length(CurrValidRecords))))
                %advance one
                lCurrRecord = CurrValidRecords(lNearInd + direction);
            else
                if direction == 1
                    lCurrRecord = CurrValidRecords(end);
                else
                    lCurrRecord = CurrValidRecords(1);
                end
            end
        end
        
        if ~isempty(lCurrRecord)
            %Call procedure to update text/slider values to this record
            SetSliderValsToRecord(0,0,lCurrRecord);
        end
    end

    function [blValue] = ToggleMenuItem(hind,event)
       %Toggles the check mark on the menu item and
       %returns the corresponding boolean value
       %If the object is a control, simply returns the value
       if strcmpi(get(hind,'type'),'uimenu')
           blValue = Enabled2Boolean(get(hind,'Checked'),1);
           set(hind,'Checked',Boolean2Enabled(blValue));
       elseif strcmpi(get(hind,'type'),'uicontrol') && strcmpi(get(hind,'Style'),'checkbox')
           blValue = get(hind,'Value');
       end
    end

    %Handlers for the Cartesian Inset plot
    % 
    function [hInsetPlotLocal,hBackBoxLocal] = LaunchInset(hWindMain)
        %Create a hidden axes and back plot
        %hWindMain = handle to the main Plot window
        %fontsizeMain = fontsize of the main window
        hInsetPlotLocal = axes('Parent',hWindMain,'box','on','visible','off');
        hBackBoxLocal = annotation(hWindMain,'rectangle','FaceColor',[0.8 0.8 0.8],'visible','off');
    end

    function ToggleInset(hInd,event)
        cVal = ToggleMenuItem(hInd);
        UpdateLegend;
        if cVal==1
            %Show the Pareto Plot
            UpdateInset;
        else
            %Hide the objective
            set([hInsetPlot hBackBox],'Visible','off');
            set(get(hInsetPlot,'Children'),'Visible','off');
        end 
    end

    function UpdateInset
       %Reposition the inset cartesian plot, update traces on it, and make
       %the ticks look nice
       
       %Check conditions are correct to proceed       
       if (nO~=2) || (strcmpi(get(uimShowInset,'Checked'),'off'))
           return;
       end
       
       hold off
       %main_fig = findobj(hWindReturn,'Type','axes');
       %ax=cell2mat(get(main_fig,'Position'));
       ax = get(plot2,'Position');
       
       inset_size = 0.2;
       InnerPos = [ax(1,1)+0.65*ax(1,3) ax(1,2)+ax(1,4)-0.65*inset_size 0.82*inset_size 0.75*inset_size]; %0.8

       set(hInsetPlot,'Position',InnerPos,'fontsize',FontSize-4);
       delete(get(hInsetPlot,'Children'));
       hold on

       cMarkers = {'.' 'o' '+' '^' 's' 's' 's' 's'};
       vSize = FontSize -12;
       %Plot each checked group in a separate color that corresponds to it's
       %renderings on the full parallel coordinate plot
       hParetoPlots = zeros(nU,1);
       for i=1:nU
           %Pull in the data for the group
           rToUse = strcmpi(vGroup,vUnique(i)); %Rows for the current group 
           %Pull in the renderings (particularly color) for the group
           [lIndToUse,lLineWidth,mColorsDecs,mColorsObjs,strVis,cMarker] = GetGroupRenderings(i,cMarkers);
           strVis = Boolean2Enabled(get(cbGroupChecks(i),'Value'));

           %Also plotted in plot units!
           hParetoPlots(i) = plot(hInsetPlot, mPlot(rToUse,2),mPlot(rToUse,1),'marker',cMarker,'MarkerSize',vSize,'color',mColorsDecs,'markerfacecolor',mColorsDecs,'linestyle','none','visible',strVis);
       end

       %Use the original axis labels, but break at the (units)
       cAxisLabels = vObjLabels;

       set(hInsetPlot,'fontsize',FontSize-6,'ycolor',[.737 0 0.737],'xcolor',[.737 0 0.737]); 
       cAxisProps = {'YLim', 'YTick', 'YTickLabel';'XLim', 'XTick', 'XTickLabel'};
       %Grab axis info from the left axis of the main plot 
       TicksAllUse = get(plot2,cAxisProps(1,:));

       %Format the tick labels on each axis
       for i=1:size(cAxisProps,1)       
           if (i==BaseAxis(1)) || ((vMults(i)==vMults(BaseAxis(1))) && all(mTransformToOrig(:,1) - mTransformToOrig(:,2) <= [1e-6;1e-6],1))
               %Use the y ticks on the main plot
               cTickCurr = TicksAllUse{2};
               cTickLabelCurr = TicksAllUse{3};
           else
               %Read in the from the ticks already on the inset plot
               cTickTemp = get(hInsetPlot,cAxisProps(i,:));
               cTickCurr = cTickTemp{2};
               cTickLabelCurr = cTickTemp{3};
           end

           %Convert from Cell to Double
           if ischar(cTickLabelCurr)
              cTickLabelCurr =  str2num(cTickLabelCurr);
           elseif iscell(cTickLabelCurr)
              cTickLabelCurr =  cellfun(@str2num,strrep(cTickLabelCurr,',',''));
           end

           %Check and if needed expand ticks to include at least 4
           if length(cTickCurr) < 4
               cTickCurr = [cTickCurr(1):(cTickCurr(end)-cTickCurr(1))/4:cTickCurr(end)];
               cTickLabelCurr = [cTickLabelCurr(1):(cTickLabelCurr(end)-cTickLabelCurr(1))/4:cTickLabelCurr(end)];
           end

           %Clip the labels to the min/max for the axis
           mLims = [min(mPlot(:,i)); max(mPlot(:,i))];
           lNumTicks = length(cTickCurr);
           iInds = [1:lNumTicks];
           vLess = iInds(((cTickCurr<mLims(1)) + (circshift(cTickCurr',lNumTicks-1)'>=mLims(1))==2)); vGreater = iInds(((cTickCurr>mLims(2)) + (circshift(cTickCurr',1)'<=mLims(2))==2));

           cLims = [cTickCurr(vLess) cTickCurr(vGreater)];
           cTickCurr = cTickCurr(vLess:vGreater);
           cTickLabelCurr = cTickLabelCurr(vLess:vGreater);

           %Check and if needed expand ticks to include at least 4
           if length(cTickCurr) < 4
               cTickCurr = [cTickCurr(1):(cTickCurr(end)-cTickCurr(1))/4:cTickCurr(end)];
               cTickLabelCurr = [cTickLabelCurr(1):(cTickLabelCurr(end)-cTickLabelCurr(1))/4:cTickLabelCurr(end)];
           end

           vMultPareto = floor((log(ConvertPlot2DataUnits(cTickCurr(end),mTransformToOrig(:,i)))-log(cTickLabelCurr(end)))/log(10));

           sTickLabel = ThousandSep(cTickLabelCurr);

           %Update the axis
           set(hInsetPlot,cAxisProps{i,1},cLims,cAxisProps{i,2},cTickCurr,cAxisProps{i,3},sTickLabel);  

           if (~isempty(strfind(cAxisLabels{i},'('))) || (vMultPareto~=0)
               [cFir cSec] = strsplit(cAxisLabels{i},'(');
               if length(cFir)==1
                   cFir{2} = ')';
               end
               if i==1 %add newline between label and units for y axis
                   cAxisLabels{i} = sprintf('%s\n(10^{%d} %s',cFir{1},vMultPareto,cFir{2});
               else
                   cAxisLabels{i} = sprintf('%s (10^{%d} %s',cFir{1},vMultPareto,cFir{2});
               end
           end
       end

       ylabel(cAxisLabels{1},'FontSize',FontSize-2,'color',[.737 0 0.737]);
       xlabel(cAxisLabels{2},'FontSize',FontSize-2,'color',[.737 0 0.737]);

       %Create a blank box to go behidnd the incoming figure
       set(hInsetPlot,'visible','on');
       vOuterPos = get(hInsetPlot,'OuterPosition');    
       %hBackBox = annotation(hWindReturn,'rectangle',vOuterPos,'FaceColor',[0.8 0.8 0.8]);
       set(hBackBox,'position',vOuterPos,'Visible','on');
 
       drawnow; pause(0.1); %let the gui catch up
       uistack(hInsetPlot,'top');
    end

%%  Generate Solution Input Box Callback Functions

    function SetFromGenBoxes(hind,event)
       % Call back that enables/disables #Samples, Gams File controls based
       % on the settings in the Generate Type,Using dropdown controls
       
       % Selecting Random Sample forces Matlab setting on Using checkbox
       % and makes the # Samples box active.
       % Selecting Gams enables the Gams File controls
       
       % Determine the settings       
       lValType = get(cbGenerateType, 'value');
       sEnableSample = 'off';
       sEnableGams = 'off';
                   
       if lValType == 3
           sEnableSample = 'on';
           set(cbGenerateMethod,'Value',2);
       end
       
       if lValType == 4
           set(cbGenerateMethod,'Value',3);
       end
       
       lValUse = get(cbGenerateMethod,'value');
       
       if lValUse == 3
           sEnableGams = 'on';
       end

       %Implement the settings      
       set(lblGamsFile,'Enable',sEnableGams);
       set(txtGamsFile,'Enable',sEnableGams);
       set(lblNumSamples,'Enable',sEnableSample);
       set(txtNumSamples,'Enable',sEnableSample);             
    end

    if hWind ~= hWindReturn
        %Trigger a resize so the plot2 and ax2 axes sync correctly
        set(hWindReturn,'Units','normalized');
        set(hWindReturn,'OuterPosition',[0.02 0.04 .96 .92]);
        set(hWindReturn,'Units','characters');
    end

    hold off
end

%% Callback Functions for Controls on Plot that generally require calling the function again after running

function PrintToPDF(hind,event,hWindCurr,mData,nO)
  %Prints the current figure window to a PDF file specified by the user.
  %Needs to delete several controls on the figure so the PDF prints
  %correctly. Then calls the function again
  
  % INPUTS
  % hWindCurr = handle to the figure
  
  %Grab the current figure settings
  [vParams,hControls] = ReadControls('PrintToPDF',hWindCurr);
  [lShowControls,lStartTab] = aGetField(vParams,{'ShowControls','StartTab'});
  [mObjs,mDecs] = SplitMatrix(mData,nO);
  %Prompt the user for the file name
  sFileName = inputdlg('Name of PDF file:','Print to PDF',1,{'output.pdf'});
  if isempty(sFileName)
     warning('No file name specified. Exiting')
     return
  end
  
  %Disable the mouse callbacks
  set(hWindCurr,'WindowButtonMotionFcn', {});
  set(hWindCurr,'WindowButtonDownFcn', {});
  
  %We need to delete several of the panel controls so the figure prints correctly
  hPanels = findall(hWindCurr,'Type','uipanel');
  
  if lShowControls==1
    %Only delete the hidden panels
    hPanels = hPanels(strcmpi(get(hPanels,'Visible'),'off')==1);
  end
  delete(hPanels);

  %Also selete the box behind the pareto plot if it exists
   hBoxes = findall(hWindCurr,'Type','hggroup');
   if ~isempty(hBoxes)
        delete(hBoxes)
   end 
 
  %Print to pdf
  PrintNearOptFigToPdf(hWindCurr,sFileName{:})
  %Recall the plotter since we deleted controls
  nearoptplotmo2(mObjs,mDecs,vParams{:});  
end

function SetVisibilityAll(cPanel,Visibility)
   % Sets visibility on the all the children controls on cPanel as well as
   % all the children panels and the children controls on those panels.
   % Recurrsively calls the function until no more panels are reached
   %  cPanel - handle of start panel
   %  Visibility - either "on" or "off" to indicate the visibility
   
    hChildren = get(cPanel,'children');
    nC = max(size(hChildren));
    
    for i=1:nC
        switch(get(hChildren(i),'Type'))
            case 'uicontrol'
                %Only change visibility for controls that are vPos > 0
                vPosTemp = get(hChildren(i),'Position');
                if vPosTemp(2)>0
                    set(hChildren(i),'Visible',Visibility);
                end
            case 'uipanel'
               SetVisibilityAll(hChildren(i),Visibility);
        end
    end
end

function [varargout] = aGetField(inArray,fields)
%   searches the cell array inArray (a 1 x D vector) for field 'field' then
%       returns the value in the field following the field name. Assumes inArray is orderd:
%       {field1name field1value field2name field2value ....}
%       returns as many values as in the cell array fields

%   returns empty if no field or field name is in last entry of the cell
%   array
    if ischar(fields)
        fields = {fields};
    end
    
    nF = length(fields);
    varargout = cell(1,nF);
    
    %fields
    %inArray

    for i=1:nF
        varargout{i} = inArray{circshift(strcmpi(inArray,fields{i})',1)};
    end
end

function [OutArray] = aSetField(inArray,varargin)
%   varargin is a list of field names followed by values, e.g.,
%   field1,field1value,field2,field2value
%
%   Inserts the values into the inArray in the field following the name.
%   If the field does not exist, appends the field name and value to the
%   end

    nF = length(varargin);
    nFInRay = length(inArray);
    OutArray=inArray;
    lNewFields=0;
    
    for i=1:2:nF-1
        %find field in the OutArray
        pos = circshift(strcmpi(inArray,varargin{i})',1);
        
        if sum(pos)==0
            %field doesnot exist append to end
            lNewFields = lNewFields+1;
            OutArray{nFInRay+2*lNewFields-1} = varargin{i};
            OutArray{nFInRay+2*lNewFields} = varargin{i+1};
        else
            %insert into existing field
            OutArray{pos}=varargin{i+1};
        end
    end 
end

function [DataVal]= ConvertPlot2DataUnits(PlotVal,ConversionFactors,blData2Plot)
  % Converts a PlotVal to a DataVal using the linear coefficients in
  % Conversion Factors. Row 1 is the slope, Row 2 the intercept.
  % Optional blData2Plot does the reverse transformation (i.e., the user
  % actually supplied a DataVal and wants a PlotVal in return)

    if nargin<=2
        blData2Plot = 0;
    end
    if blData2Plot
        DataVal = (PlotVal-ConversionFactors(2))/ConversionFactors(1);
    else
        DataVal = ConversionFactors(1)*PlotVal+ConversionFactors(2);
    end
end

function [A,B] = SplitMatrix(C,nCol)
    %Splits matrix C at the nCol specified column so that columns 1..nCol
    %are in Matrix A and nCol+1..n are in Matrix B
    
    [m,n]=size(C);
    
    if nCol>=n
        A = C;
    elseif nCol==0
        A=zeros(m,0);
        B=C;
    else
        A = C(:,1:nCol);
        B = C(:,nCol+1:n);
    end
end

function [mMatrixReturn] = CompactToFull(mMatrix,vChecks,vCheckedValues,ReverseMap)
    %Takes the m x nS matrix mMatrix and returns a fuller m x n matrix with
    %the values for the n - nS columns not in the input matrix filled from
    %the vector vCheckedValues. vChecks is a logical boolean vector (1 x
    %n). Values of indicate draw column values from the vector
    %vCheckedValues. ReverseMap = 1 specifies the reverse transformation (Full to Compact)
    
    %INPUTS
    % mMatrix = m x nS matrix of input data
    % vChecks = 1 x n vector of booleans. Value of 1 indicates the column
    %   is fixed and column values should be filled from the ith element in
    %   vCheckedValues
    % vCheckedValues = 1 x n vectors of values used to fill checked columns
    % ReverseMap = boolean, takes value of one to indicate reverse mapping
    % (Full to Compact). 0 or omitted indicates Compact To Full.
    
    %OUTPUTS
    % mMatrixReturn = m x n matrix comprised of the input data and values
    %      filled from vCheckedValues. If ReverseMap=1, then mMatrixReturn
    %      is a m x nS matrix where nS = the number of unchecked (0-value)
    %      entries in vChecks and valus are from those columns of the 
    
    %Handle optional parameter
    if nargin < 4
        ReverseMap = 0;
    end
    
    [m nS] = size(mMatrix);
    nC = length(vChecks);
    nC2 = length(vCheckedValues);
    
    if ReverseMap == 1
        %Full to compact
        if nS~=nC
            error(sprintf('nearoptplot2:CompactToFull-mMatrix and vChecks need to have the same number of columns. Currently: %d and %d',nS,nC));
            return
        end
        mMatrixReturn = mMatrix(:,vChecks==0);
    else
        %Compact to full
        if (nS+sum(vChecks==1)~=nC) || (nC~=nC2)
            error(sprintf('nearoptplot2:CompactToFull-Number of columns in mMatrix plus # of checked columns must equal size of vChecks. Also, vCheckedValues and vChecks must have same size. Currently, mMatrix%d, Checked: %d, vChecks: %d, vCheckedVals: %d',nS,sum(vChecks==1),nC,nC2));
            return
        end
        mMatrixReturn = zeros(m,nC);
        mMatrixReturn(:,vChecks==0) = mMatrix(:,:);
        mMatrixReturn(:,vChecks==1) = repmat(vCheckedValues(vChecks==1),m,1);
    end
end

function [AllowableDeviationVal, AllowableDeviationValPU] = GetCheckAllowableDeviation(txtAllowableDeviation,strCallFunction,TickSize)
    %returns the Allowable deviation in the checkbox. Checks that the
    %value is numerical and returns the absolute value
    %If the optional parameter TickSize is passed, additionally calculates
    %and returns AllowableDeviationVal in Plot Units.
    
    AllowableDeviationVal = get(txtAllowableDeviation,'String');

    %error checking on AllowableDeviationValue
    if ~isnumeric(str2num((AllowableDeviationVal)))
            warning(['nearoptmo2: ',strCallFunction],'AllowableDeviation is not a numerical value. Defaulting to zero to continue')
            AllowableDeviationVal = 0;
    else
            AllowableDeviationVal= abs(str2num(AllowableDeviationVal));
    end
    
    if nargin>2
        AllowableDeviationValPU=AllowableDeviationVal*TickSize;
    else
        AllowableDeviationValPU = NaN;
    end
end

function [vSliderValuesDU vSliderSteps vSliderValuesPU vMins vMaxs] = ReadSliderValues(sSliders,vFixedValues,mConvert)
    %returns the vector of sSliders values. Where a sSlider is not defined,
    %instead substitutes the predefined fixed value.
    %Note vFixedValues are provided in Data units, while sliders are in
    %plot units.
    %Outputs slider values in both data units (vSliderValuesDU) and plot
    %units (vSliderValuesPU)
    %Use the coefficients mConvert to convert between the two.
    %Also outputs the slider min and max values in plot units (vMins and
    %   vMaxs)
    
    vSliderValuesDU = vFixedValues;
    vSliderValuesPU = vSliderValuesDU;
    vSliderValuesPU(:)=0;
    vMins = vSliderValuesPU;
    vMaxs = vMins;
    vSliderSteps = vFixedValues;
    n = max(size(vSliderValuesDU));
    
    for i=1:n
        if sSliders(i)>0
            vSliderValuesPU(i) = get(sSliders(i),'Value');
            sSlideSteps = get(sSliders(i),'SliderStep');
            vMins(i) = get(sSliders(i),'Min');
            vMaxs(i) = get(sSliders(i),'Max');
            vSliderSteps(i) = sSlideSteps(1);
            vSliderValuesDU(i)=ConvertPlot2DataUnits(vSliderValuesPU(i),mConvert(:,i));
        else
            vSliderValuesPU(i)=ConvertPlot2DataUnits(vSliderValuesDU(i),mConvert(:,i),1);
        end
    end
    %fprintf('Inside Read Slider Values\n');
    %vSliderValues
end

function [vCheckValues] = ReadCheckboxValues(cbChecks)
    % returns the vector of cbChecks box values into the vector vCheckValues
    [m n] = size(cbChecks);
    vValues = zeros(1,n);
    
    for i=1:n
        vValues(i)=get(cbChecks(i),'Value');
    end    

    vCheckValues = vValues;
end

function [outargs,hControls,mObjs,mDecs] = ReadControls(CallFunction,hWindCurr)
    %Queries select entries on hControls and updates corresponding values
    %in varargin and returns as cell array outargs.
    %Also returns:
    %   hControls = cell array list of {label, handle, label, handle, ...}
    %       pairs to the controls
    %   mObjs = m x nO matrix of stored objective function values (nO=number of
    %       objectives)
    %   mDecs = m x nD matrix of stored decision variable values (nD =
    %       number of decision variables)
    %
    % For the case where no~=size of tolerances (added new objective(s), controls are of a different size, instead update values from
    % the variable argument list
    
    % Get the handles for the controls we will need
    hControls = getappdata(hWindCurr,'hControls');
    varargin = getappdata(hWindCurr,'varargs');
    mConvert = getappdata(hWindCurr,'mConvert');
    mObjs = getappdata(hWindCurr,'mObjs');
    mDecs = getappdata(hWindCurr,'mDecs');
    
    nO = size(mObjs,2);
    
    [txtTolerances, txtNumSamples, txtObjSamples, cbChecks, txtAllowableDeviation, sSliders, txtGroupOrders,txtGroupNames,txtGroupThicks,cbGroupChecks,txtGamsFile,cTabs,txtCurrRec,cbGenType,cbGenMethod,cbGroupToHighlight,txtFontSize] = aGetField(hControls,{'Tolerance','NumSamples','ObjSamples','AxisChecked','AllowableDeviation','Sliders','GroupOrders','GroupNames','GroupThicks','GroupChecks','GamsFile','Tabs','CurrentRecord','cbGenType','cbGenMethod','HighlightGroup','txtFontSize'});
    [uimHideSliders,uimHideChecks, uimHideControls, uimShowCurrRecord, uimShowFiltered, uimShowObjsDiffColor, uimGroupLabels, uimShowInset, uimShowLegendOptions] = aGetField(hControls,{'HideSliders','HideCheckboxes','HideControls','ShowCurrRecord' 'ShowFilteredAlts' 'ShowObjsDiffColor' 'ShowGroupLabels' 'ShowInsetPlot' 'LegendOptions'});
    
    %Update the plot
    drawnow; pause(0.1);
    %Group data
    mGroupData = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks,txtGroupThicks);
    %Text boxes
    %Tolerance(s)
    if nO ~= length(txtTolerances) % un-equal size because added objective(s), update values from values preloaded into varargin
        [ToleranceValues,mGroupData,vFixed]  = aGetField(varargin,{'Tolerance','mGroupData','vFixed'});
    elseif (nO == 1) % Single box
        ToleranceValues = str2num(get(txtTolerances,'String'));
    else
        %Multiple boxes for the multiple objectives
        ToleranceValues = cellfun(@str2num,(get(txtTolerances,'string')));
        ToleranceValues = ToleranceValues';
    end
    
    if any(isempty(ToleranceValues))
        error('Un-acceptable near-optimal tolerance input. Must be numeric. Try again');
        return;
    end
    
    sGamsFile = get(txtGamsFile,'String');   
    if ~strcmpi(sGamsFile,'') && (exist(sGamsFile,'file') == 0)
        %File does not exist
        error(sprintf('nearoptplotmo2 (%s) - File %s does not exist',CallFunction,sGamsFile));
    end
    
    NumSamples = GetCheckAllowableDeviation(txtNumSamples,CallFunction);
    ObjSamplesPercent = GetCheckAllowableDeviation(txtObjSamples,CallFunction);
    if ObjSamplesPercent > 100
        ObjSamplesPercent = 100;
    end
    
    NewFontSize = str2num(get(txtFontSize,'String'));
    if isempty(NewFontSize)
        NewFontSize = aGetField(varargin,{'FontSize'});
        warning(sprintf('Font Size %s is not a valid input. Reverting to prior setting.',get(txtFontSize,'String')));
    end
    
    %Start Tabs
    lTabInds = [1:length(cTabs)];
    StartTab = lTabInds(strcmpi(get(cTabs,'visible'),'on'));      %  StartTab = GetAllValues(cTabs,'visible','on');
    %Current Record
    lCurrRecord = str2num(get(txtCurrRec,'String'));
    %Generate Type and method
    GenType = get(cbGenType,'value');
    GenMethod = get(cbGenMethod,'value');
    %Hide sliders and Hide controls from menu selections
    HideSliders = Enabled2Boolean(get(uimHideSliders,'Checked'));
    HideCheckboxes = Enabled2Boolean(get(uimHideChecks,'Checked'));
    ShowControls = Enabled2Boolean(get(uimHideControls,'Checked'),1);
    ShowCurrRecord = Enabled2Boolean(get(uimShowCurrRecord,'Checked'));
    ShowFilteredAlts = Enabled2Boolean(get(uimShowFiltered,'Checked'));
    ShowObjsDiffColor = Enabled2Boolean(get(uimShowObjsDiffColor,'Checked'));
    ShowGroupLabels = Enabled2Boolean(get(uimGroupLabels,'Checked'));
    ShowInsetPlot = Enabled2Boolean(get(uimShowInset,'Checked'));
    
    LegendOptions = get(uimShowLegendOptions,'Label');
    LegendChecks = get(uimShowLegendOptions,'checked');
    ShowLegend = LegendOptions(strcmpi(LegendChecks,'on'));
    
    %Group to highlight
    lGroupToHighlight = get(cbGroupToHighlight,'Value');
    if lGroupToHighlight <= 1
        GroupToHighlight = '0';
    else
        cGroupNames = get(cbGroupToHighlight,'String');
        GroupToHighlight = deblank(cGroupNames(lGroupToHighlight,:));
    end
        
    %Read in the variable arguments
    [vFixedVals TickSize] = aGetField(varargin,{'vFixedVals' 'TickSize'});   
    
    [AllowableDeviationValue, AllowableDeviationValPU] = GetCheckAllowableDeviation(txtAllowableDeviation,CallFunction,TickSize);
    if nO==length(txtTolerances) %update vFixedVals from control settings
        vFixed = ReadCheckboxValues(cbChecks);   %Axes check boxes
        [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    end
    
    %Update the values in the field list
    outargs = aSetField(varargin,'Tolerance',ToleranceValues,'AllowableDeviation',AllowableDeviationValue,'AllowableDeviationPU',AllowableDeviationValPU, ...
            'NumSamples',NumSamples,'ObjSamplesPercent',ObjSamplesPercent,'vFixed',vFixed,'vFixedVals',vFixedVals,'mGroupData',mGroupData,'sGamsFile',sGamsFile, ...
            'mConvert',mConvert,'StartTab',StartTab,'CurrentRecord',lCurrRecord,'GenerateType',GenType,'GenerateMethod',GenMethod, ...
            'HideSliders',HideSliders,'HideCheckboxes',HideCheckboxes,'GroupToHighlight',GroupToHighlight,'ShowControls',ShowControls,'FontSize',NewFontSize, ...
            'nO',nO,'ShowCurrRecord',ShowCurrRecord, 'ShowFilteredAlts',ShowFilteredAlts, 'ShowObjsDiffColor',ShowObjsDiffColor, ...
            'ShowGroupLabels', ShowGroupLabels, 'ShowInsetPlot',ShowInsetPlot,'ShowLegend',ShowLegend);
end

function [ReturnVal] = GetCheckTxtValue(strCallFunction, txtCntl)
%error checking on the input to txtError to make sure it's numeric

   TempVal = get(txtCntl,'string');
    
    if ~isnumeric(str2num((TempVal)))
        warning(['nearoptmo2: ', strCallFunction],['Text input for ', get(txtCntl,'name'), ' is not a numerical value. Defaulting to zero to continue'])
        ReturnVal = 0;
    else
        ReturnVal = abs(str2num(TempVal));
    end
end

function [gNewOrder] = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks,txtGroupThicks)
    % Takes the input on txtGroupOrders and txtGroupNames boxes all the group info (txtGroupNames, cbGroupChecks, txtGroupThicks)
    % Returns a cell matrix of the newly sorted Group Names and Checked
    % values.

    %Read in the new order
    lNumGroups = length(txtGroupOrders);
    gOrder = cell(lNumGroups,4);
    for i=1:lNumGroups
        gOrder(i,1) = {get(txtGroupOrders(i),'string')};
        gOrder(i,2) = get(txtGroupNames(i),'string');
        gOrder{i,3} = get(cbGroupChecks(i),'value');
        gOrder{i,4} = str2num(get(txtGroupThicks(i),'string'));
    end
    
    %gOrder;

    [gSort, gIC] = sortrows(gOrder,1);
    
    %gSort;

    %Reassign the sorted variables
    gNewOrder = gSort(:,[2:4]);
end

function HighlightGroup(hObj,event,mResults,nO,hWindCurr)
    %Assigns the selected group in the HighlightGroup popup control as the
    %new highlighted groups. Moves this group to the final position so it
    %plots on top.

    [vParams,hControls] = ReadControls('ShowGroups',hWindCurr);
    [mGroupInfo,GroupToHighlight] = aGetField(vParams,{'mGroupData','GroupToHighlight'});
    [txtOrders,cbGroupChecks] = aGetField(hControls,{'GroupOrders','GroupChecks'});
    
    if ~strcmpi(GroupToHighlight,'0')
        %A group is selected to highlight. Set this group order to be
        %highest so it plots last and on top. Check the group so it is plotted.
        lGroupHigh = find(strcmpi(mGroupInfo(:,1),GroupToHighlight),1,'first');
        set(txtOrders(lGroupHigh),'String',num2str(size(mGroupInfo,1)+1));
        set(cbGroupChecks(lGroupHigh),'Value',1);
    end
    
    ReorderGroups(0,0,mResults,nO,hWindCurr);
end


function RemoveGroups(hObj,event,mResults,nO,hWindCurr,lAction)
    %Reassigns checked groups based on the lAction selected
    % 1 - re-assigns all checked groups to the first checked group
    % 2 - deletes the checked groups (permanently removes)
    
    %varargin_start = getappdata(hWindCurr,'varargs');
    [varargin_fresh,hControls] = ReadControls('RemoveGroups',hWindCurr);
    [vGroup iOpt mGroupData OptSolRow] =  aGetField(varargin_fresh,{'vGroup' 'GroupToHighlight' 'mGroupData' 'OptSolRow'});

    [cbGroupChecks] = aGetField(hControls,{'GroupChecks'});
  
    vChecked = cell2mat(mGroupData(:,2));
    vUniques = mGroupData(:,1);

    if sum(vChecked) == 0
        warning('Need to check at least one group before proceeding.')
        return
    end
    
    m = size(vGroup,1);
    nG = length(vUniques);
    vGroupNew = vGroup;
    
    switch lAction
        
        case 1 %re-assign all checked groups to the first checked group
            
            if sum(vChecked) < 2
                msgbox('Need to check at least two groups to merge.','Warning')
                return
            end
    
            %Grab the index of the first checked group
            [gTemp, fInd] = GetFirst(cbGroupChecks,vUniques,1);
            gInds = [1:nG];
            gIndKeep = gInds((vChecked'==0) + (gInds==fInd) >= 1);
            gIndsMerge = gInds((vChecked'==1) + (gInds>fInd) == 2);
            gRevertTo = vUniques(fInd);

            for i=gIndsMerge
                vCurrG = strcmpi(vGroupNew,vUniques(i)); %indices of current group
                vGroupNew(vCurrG==1) = repmat(gRevertTo,size(sum(vCurrG),1),1); 
            end
            
            mGroupDataNew = mGroupData(gIndKeep,:);
            mResultsNew = mResults;
            OptSolRowNew = OptSolRow;

        case 2 %delete the checked groups
            %Remove those rows from vGroup and mResults
            mDeleteGroups = mGroupData(vChecked==1,1)';
            sAdd = '';
            if length(mDeleteGroups)>1
                sAdd = 's';
            end
            
            sGroupList = strjoin(mDeleteGroups,', ');
            sResp = questdlg(['Delete data associated with the group',sAdd,': ', sGroupList,'? (You can not later retrieve in the plot)'], 'Delete groups?','Yes','No','Yes');
            drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
            if strcmpi(sResp,'No')
                return;
            end
            
            mGroupDataNew = mGroupData(vChecked==0,:);
            RowsToRetain = zeros(m,1);
            if isnumeric(vGroup)
                for i=1:m
                    RowsToRetain(i) = any(vGroup(i)==mGroupDataNew(:,1));
                end
            else
                for i=1:m
                    RowsToRetain(i) = any(strcmpi(vGroup(i),mGroupDataNew(:,1)));
                end
            end
            
            vGroupNew = vGroup(RowsToRetain==1);
            mResultsNew = mResults(RowsToRetain==1,:);
            %Update the optimal solution row
            if all(RowsToRetain(OptSolRow)) && (~isempty(OptSolRow))
                %Optimal solution is retained. Determine new row #
                %in retained rows. This is simply the count of the retained rows
                %up to current OptSolRow
                for i=1:nO
                    OptSolRowNew(i) = sum(RowsToRetain(1:OptSolRow(i),1));
                end
            else  %Optimal solution is not retained; set back to default
                warning('You deleted a group that had an optimal solution as a member. Will not be able to recalculate near-optimal constraint')
                OptSolRowNew = [];
            end      
    end 
    
    %Revert to default of show all retained groups
    for i=1:size(mGroupDataNew,1)
        mGroupDataNew{i,2} = 1;
    end
            
    varargout = aSetField(varargin_fresh,'vGroup',vGroupNew,'mGroupData',mGroupDataNew,'OptSolRow',OptSolRowNew);
    [mObjs, mDecs]=SplitMatrix(mResultsNew,nO);
    nearoptplotmo2(mObjs, mDecs,varargout{:}); 
end

function ReorderGroups(hObj,event,mResults,nO,hWindCurr)
    %Reorders the groups by the entries in txtGroupOrders
    %Simply query the controls and re-call the plot function    
    [vararg_out,hControls] = ReadControls('ReorderGroups',hWindCurr);
    [mObjs, mDecs] = SplitMatrix(mResults,nO);   
    nearoptplotmo2(mObjs, mDecs,vararg_out{:});     
end

function [sNewGroupName, mGroupDataNewSort] = GenerateNewGroup(hWindCurr,dNewGroup,sAppend)
    %Generates a new group name based on the entries in the
    %GroupOrders,GroupNames, and GroupChecks controls
    % This name is a derivative of the first checked group or the dNewGroup value passed
    
    %INPUTS
    % hWindCurr = handle of the originating figure
    % dNewGroup = optional parameter when included is the new group name
    %               given
    % sNamePart = optional parameter when included is a fixed part appended
    %    as part of the new name
    %    
    %OUTPUTS
    % sNewGroupName = the new group name generated
    % mGroupDataNewSort = matrix of the new group names and data in the
    %                       order they should appear in the list
    
    [vararg_out,hControls] = ReadControls('GenerateNewGroup',hWindCurr);
    [mGroupDataSort, vGroup] = aGetField(vararg_out,{'mGroupData' 'vGroup'});
    cbGroupChecks = aGetField(hControls,{'GroupChecks'});
    
    %Work on the grouping. The new group will appear with a child name just below the
    %first checked group in the resorted list
    
    %Find the first checked group
    [gValue, gDex] = GetFirst(cbGroupChecks,mGroupDataSort(:,1),1);

    if iscell(vGroup)          
        if nargin == 2
           % Add the specified name as a new full name below the first checked group
           sNewGroupName = dNewGroup;
        else
           %auto generate the name from the 1st checked name. Append the name part
           sNewGroupName = sprintf('%s_%s',mGroupDataSort{gDex,1},sAppend);
        end                 
       
       gSortVal = gDex + 0.1;
       sNewGroupNameBase = sNewGroupName;
       
        i = 0;
        while (sum(strcmpi(sNewGroupName,mGroupDataSort(:,1))) > 0) && (i<10)
           %The new name already exists, add a number increment to the end 
           i=i+1;
           sNewGroupName = sprintf('%s%d',sNewGroupNameBase,i+1);
        end

        gSortVal = gDex+i+0.1;
    else
        %auto generate the group name as a numeric
        sNewGroupName = gValue+1/10;

        i=0;
        while (sum(sNewGroupName==mGroupDataSort(:,1))> 0) && (i<10)
            %The index already exists, increment
            i=i+1;
            sNewGroupName = gValue+(i+1)/10;
        end
        gSortVal = sNewGroupName;
    end   

    %Resort the GroupData with the new group below the first checked
    %mGroupDataSort
    gSortStart = cell(size(mGroupDataSort,1),1);
    for i=1:size(mGroupDataSort,1)
        gSortStart{i} = i;
    end

    mGroupDataNew = [mGroupDataSort(:,1:3) gSortStart(:); {sNewGroupName} {1} mGroupDataSort(gDex,3) {gSortVal} ]; %Line thickness stays the same as First checked group
    mGroupDataNewSort = sortrows(mGroupDataNew,4);
    
end
                                    
%function Reorder(hObj,event,mMatrix, nObjs, txtTolerance, cbChecks, sSliders,mConvert,cbReorder, cbByCat, varargin) %#ok<INUSL>
function Reorder(hObj,event,mMatrix, nObjs,hWindCurr)
    % Reorder the decision variable columns of mMatrix according to the dynamic ranges of the decision variables
    % If cbReorder is checked then put the constant non-zero decision variables first, dyanmic range variables second, and zero-valued variables third.
    % If cbReorder is NOT checked: put the variables with dynamic range
    % first, constant non-zero decision variables second, and zero-valued decision variables third
    
    % If cbByCat is checked then do the organization by the overarching
    % category type in vActCats (i.e., all 1's first, 2's second, etc...)
    
    % assumes first nObjs columns of mMatrix are objective function values
    % and are not part of the reordering
    
    % Also reorders the corresponding Fixed Axes in cbChecks, labels and
    % abreviations in vXLables and vXLabelsShort, and vActCats
       
    %[max(mMatrix(:,2:end),2); min(mMatrix(:,2:end),2)]
    
    [m,n] = size(mMatrix);
    
    [varargin,hControls] = ReadControls('Reorder',hWindCurr);
    %Read in the varargin parameters
    [vFixedVal vStep vXLabels vXLabelsShort mActCat vObjLabels] = aGetField(varargin,{'vFixedVals' 'vStep','vXLabels','vXLabelsShort','mActCat','vObjLabels'});
    [vFixed Tolerance] = aGetField(varargin,{'vFixed' 'tolerance'}); 
    
    [mS,nS] = size(mActCat);
    %Read in the check box values
    %vFixed = ReadCheckboxValues(cbChecks);
    %[vFixedVal vSliderSteps vFixedValPU] = ReadSliderValues(sSliders,vFixedVal,mConvert);
    %Tolerance = get(txtTolerance,'String');
    [cbReorder, cbByCat] = aGetField(hControls,{'Reorder' 'ByCat'});
    iOrderType = get(cbReorder,'Value');
    iUseCat = get(cbByCat,'Value');
    
    mMatrixN = mMatrix;
    vFixedN = vFixed;
    vFixedValN = vFixedVal;
    vStepN = vStep;
    vXLabelsN = vXLabels;
    vXLabelsShortN = vXLabelsShort;
    mActCatSort = mActCat;
    
    if iUseCat==1
        %reorganize columns by category
        [mActCatSort iCatSort] = sortrows(mActCat,[1:nS]);
        [cats,iCatStart,iCat] = unique(mActCatSort(:,1));
        %count columns in each category
        CatCounts = histc(iCat,unique(iCat));
        nCats = length(cats);
       
        for i=1+nObjs:n
            mMatrixN(:,i) = mMatrix(:,iCatSort(i-nObjs)+nObjs);
            vFixedN(i) = vFixed(iCatSort(i-nObjs)+nObjs);
            vFixedValN(i) = vFixedVal(iCatSort(i-nObjs)+nObjs);
            vStepN(i) = vStep(iCatSort(i-nObjs));
            vXLabelsN{i-nObjs} = vXLabels{iCatSort(i-nObjs)};
            vXLabelsShortN{i-nObjs} = vXLabelsShort{iCatSort(i-nObjs)};
        end
    else
        %there is only one category
        nCats = 1;
        CatCounts = [n-nObjs];
    end
    
    begCol = nObjs+1;
    %now recast the matrix
    mResult = mMatrixN;
    vNewXLabels = vXLabelsN;
    vNewFixed = vFixedN;
    vNewXLabelsShort = vXLabelsShortN;
    vNewFixedVal = vFixedValN;
    vNewStep = vStepN;
    mActCatNew = mActCatSort;
    
    
    for k=1:nCats  %loop through the categories  
        
        %initialize the variables that will hold column numbers for each
        %category
        vDynamic = zeros(0);
        vNonZero = zeros(0);
        vZero = zeros(0);
        
        

        %loop through the columns and classify each column
        for i=begCol:(begCol+CatCounts(k)-1) %(nObjs+1):n
            [k i begCol CatCounts(k)];
            if max(mMatrix(:,i)) > min(mMatrix(:,i))
               vDynamic = [vDynamic; i];
           elseif max(mMatrix(:,i)) > 0 
               vNonZero = [vNonZero; i];
           else
                   vZero = [vZero; i];
           end
        end

        vNonZero;
        vDynamic;
        vZero;
        
        if iOrderType == 0 %dynamic ranged variables go first
            vCurrCols = [vDynamic; vNonZero; vZero];
        else
            %constant nonzeros go first
            vCurrCols = [vNonZero; vDynamic; vZero];
        end
        
        if k==1
            vNewCols = vCurrCols;
        else
            vNewCols = [vNewCols; vCurrCols];
        end
        begCol = begCol+CatCounts(k);
    end

    %copy over the columns
    [vDm, vDn] = size(vNewCols);
    
    for i=1:vDm

       mResult(:,i+nObjs) = mMatrixN(:,vNewCols(i));
       vNewFixed(i+nObjs) = vFixedN(vNewCols(i));
       vNewFixedVal(i+nObjs) = vFixedValN(vNewCols(i));
       mNewActCat(i,:) = mActCatSort(vNewCols(i)-nObjs,:);
       vNewStep(i) = vStepN(vNewCols(i));
       vNewXLabels{i} = vXLabelsN{vNewCols(i)-nObjs};
       vNewXLabelsShort{i} = vXLabelsShortN{vNewCols(i)-nObjs};
    end

    %[max(mMatrix(:,2:end)); min(mMatrix(:,2:end)); vNewCols']
    
    size(mResult);
    mNewActCat;
    vNewXLabels;

    [mObjs, mDecs]=SplitMatrix(mResult,nObjs);
    
    %Set parameters that were manipulated to pass via varargin
    varargout = aSetField(varargin,'vFixed',vNewFixed,'vFixedVals',vNewFixedVal,'vStep',vNewStep,'vXLabels',vNewXLabels, ...
                'vXLabelsShort',vNewXLabelsShort,'mActCat',mNewActCat);
      
    nearoptplotmo2(mObjs, mDecs,varargout{:});
end

function RearrangeAxes(hObj,event,FarDir,Direction,hWindCurr,mMatrix,nObjs) %#ok<INUSL>

    % Called to rearrange or alter the order, position, and/or number of axis plotted. This can include
    % A) Prune the number of dimensions (remove checked axes)
    % B) Prune axes whose values are all zero
    % C) Move checked axes one position to the right or left
    % D) Move checked axes to the far right or far left of the decision or
    %       objective axes groups
    % 
    % INPUTS
    % Direction:
    %   +1 = move checked axes to the right
    %    0 = prune (remove) checked axes
    %   -1 = move checked axes to the left
    %
    % FarDir
    %    1 = n move to the far-most axis postion in Direction
    %    0 = move only one-axis position in Direction
    % Reads the controls cbChecked and cbZeros in hWindCurr to further
    % specify action
    %    cbChecked - checked (value of 1) and Direction of 0 == Prune/take out all the checked axes (columns in the matrix mMatrix and associated data).
    %    cbZeros is checked and Direction of  0 = Prune/take out all zero-valued axes
    
    % Repositioning is only within the objective function or decision
    % variable axes -- i.e. can't move an objective function axes into and
    % among decision variable axes
    
    % nObjs = number of objective functions in mMatrix (first nObjs columns of mMatrix)
    
    % Also reorders/prunes the corresponding Fixed Axes in cbChecks, vObjLabels, and labels and
    % abreviations in vXLables and vXLabelsShort and all other columnar
    % data passed along via hWindCurr to facillitate replotting the chart
    
    % Step 1. Figure out new order of columns (from Direction, FarDir, cbChecked, cbZeros)
    % Step 2. Reorder all the columnar data (mMatrix, vFixedVals, vStep, vXLabels, vXLabelsShort, mActCat, vObjLabels, mLims, ProbForm, cFunc
    % Step 3. Call the plot function with the reordered columnar data
    
    [m,n] = size(mMatrix);
    %calculate the number of decision variables (non objectives)
    nD = n-nObjs;
    
    [varargin_rrax,hControls] = ReadControls('RearrangeAxes',hWindCurr);
    %Read in the varargin parameters
    [vFixedVal vStep vXLabels vXLabelsShort mActCat vObjLabels] = aGetField(varargin_rrax,{'vFixedVals' 'vStep','vXLabels','vXLabelsShort','mActCat','vObjLabels'});
    [vFixed Tolerance] = aGetField(varargin_rrax,{'vFixed' 'tolerance'});
    
    [cbChecked,cbZero] = aGetField(hControls,{'PruneChecked','PruneZeros'});
    iChecked = get(cbChecked,'Value');
    iZero = get(cbZero,'Value');
    
    AxInds = [1:n];  
    DecInds = [1:nD];
    if nObjs>0
        ObjInds = [1:nObjs];
    end
    
    if Direction==0 %Prune checked/zeroed axes
        %first identify which columns to {Prune (==1) and retain (==0)
        %start by assuming we keep them all (==0)
        vPruneCheck = zeros(1,n);
        vPruneZero = zeros(1,n);

        if iChecked>0
            vPruneCheck = vFixed;
            sError = 'No axes checked. Need to check at least one axis to indicate the axes to delete.';
        elseif iZero>0
            vMax = max(mMatrix,2);
            vMin = min(mMatrix,2);

            vPruneZero = (vMax == 0) .* (vMin==0);
            sError = 'No axes have values that are all zero. No axes to delete.';
        else
            sError = 'Indicate criteria to use to delete axes (checked or have value of zero)';
        end

        vPrune = max([vPruneCheck; vPruneZero]);  
        
        if sum(vPrune)==0
            msgbox(sError);
            set(cbZero,'Value',0);
            return
        end
            
        vKeep = AxInds(vPrune==0);
        vKeepDec=DecInds(vPrune(nObjs+1:n)==0);
        sObjs = '';
        if nObjs>0
            vKeepObj=ObjInds(vPrune(1:nObjs)==0);
            lPruneObj = sum(vPrune(1:nObjs));
            lNewObjs = nObjs-lPruneObj;
            mDelObjs = vObjLabels(vPrune(1:nObjs)==1);           
            if lPruneObj>0
                sObjs = strjoin(mDelObjs,', ');
            end
        else
            vKeepObj=[];
            lPruneObj=0;
        end
        
        vRemove = AxInds(vPrune==1);
        lPruneDec = sum(vPrune(nObjs+1:n));
        
        sDecs = '';
        if (lPruneDec > 0)
            mDelDecs = vXLabels(vPrune(nObjs+1:n)==1);
            sDecs = strjoin(mDelDecs,', ');
        end
        
        %Check with user to confirm deleting     
        if lPruneObj + lPruneDec > 1
            sAdd = 'axes';
        else
            sAdd = 'axis';
        end
        
        sJoin = '';
        if (lPruneObj>0) && (lPruneDec > 0)
            sJoin = ', ';
        end

        sResp = questdlg(['Delete the following ',sAdd,': ', sObjs,sJoin,sDecs,'? (You can not later retrieve in the plot)'], 'Delete axes?','Yes','No','Yes');
        drawnow; pause(0.05);  % this innocent line prevents the Matlab hang
        if strcmpi(sResp,'No')
            return;
        end      
                
        lPrune = sum(vPrune);
        lNewDec = nD - sum(vPrune(nObjs+1:n));
        
    else
        %Move checked axes left or right
            %calcualte the destination column
        vCheckDec =  DecInds(vFixed(nObjs+1:n)==1);
        vNoCheckDec = DecInds(vFixed(nObjs+1:n)==0);
        if nObjs>0
            vCheckObj = ObjInds(vFixed(1:nObjs)==1);
            vNoCheckObj = ObjInds(vFixed(1:nObjs)==0);
        end
        
        if FarDir==1           
            %goto the extremes
            if Direction==-1 %left means checked axes come first
                vKeepDec = [vCheckDec vNoCheckDec];
                 if nObjs>0
                    vKeepObj = [vCheckObj vNoCheckObj];
                 end
            else %right means checked axes come last
                vKeepDec = [vNoCheckDec vCheckDec];
                 if nObjs>0
                    vKeepObj = [vNoCheckObj vCheckObj]; 
                 end
            end
        else
            %move only one away in direction     
            [vDummy, vKeepDec] = ShiftNonZero(vFixed(DecInds+nObjs),Direction);
            if nObjs>0
                [vDummy, vKeepObj] = ShiftNonZero(vFixed(ObjInds),Direction);
            end
            %add back in indexes for columns that were not selected
        end
        
        if nObjs>0
            vKeep = [vKeepObj vKeepDec+nObjs];
        else
            vKeep = vKeepDec+nObjs;
        end
        
        lNewDec = nD;
        lNewObjs = nObjs;
    end
    
    %vKeep;
    %vKeepObj;
    %vKeepDec;
    
    %Read in the varargin parameters that will be maniputed
    [vFixedVal vStep vXLabels vXLabelsShort mActCat vObjLabels mLims ProbForm cFunc BaseAxis Tolerance AllowableDeviation OptSolRow NearOptConstraint blShowInset] = aGetField(varargin_rrax,{'vFixedVals' 'vStep','vXLabels','vXLabelsShort','mActCat','vObjLabels','mLims', 'ProbForm','cFunc','BaseAxis','Tolerance','AllowableDeviation','OptSolRow' 'NearOptConstraint' 'ShowInsetPlot'});
    
    %[vFixedVal vSliderSteps vFixedValPU] = ReadSliderValues(sSliders,vFixedVal,mConvert); 
    %AllowableDeviation = GetCheckAllowableDeviation(txtAllowableDeviation,'Prune');
 
    [mS,nS] = size(mActCat);
     
    if nObjs==0
        vLabelsAll = vXLabels;
    else
        vLabelsAll = {vObjLabels{:} vXLabels{:}};
    end
        
    %retain the columns to keep
    [mMatrixN, vFixedN, vFixedValN, vStepN,vLabelsAllN, mLimsN] = ReorderCols(vKeep,mMatrix,vFixed,vFixedVal',vStep,vLabelsAll,mLims);
    %retain objective-function specific columns
    [ToleranceN,OptSolRowN,NearOptConstraintN] = ReorderCols(vKeepObj,Tolerance,OptSolRow,NearOptConstraint);
    
    vXLabelsShortUse = vXLabelsShort;
    
    if size(vKeepDec) ~= size(vXLabelsShort)
        vXLabelsShortUse = vXLabelsShort';
    end
    
    [vXLabelsShortN,mActCatN] = ReorderCols(vKeepDec,vXLabelsShortUse,mActCat');

    mActCatN = mActCatN';
    
    if ~isempty(ProbForm)
       [AineqN, AeqN, lbN, ubN] = ReorderCols(vKeepDec,ProbForm.Aineq,ProbForm.Aeq,ProbForm.lb',ProbForm.ub');
        
       if lNewDec ~= nD
            %Update the right hand side of the inequality and equity
            %constraints
            bineqN = ProbForm.bineq - ProbForm.Aineq*(vPrune(nObjs+1:n).*vFixedVal(nObjs+1:n)')';
            beqN = ProbForm.beq -  ProbForm.Aeq*(vPrune(nObjs+1:n).*vFixedVal(nObjs+1:n)')';
        else
           bineqN = ProbForm.bineq;
           beqN = ProbForm.beq;
       end

       
       %Remove rows in the constraints that correspond to near-optimal
       %constraints of removed objectives
        if lNewObjs ~= nObjs
            lCons = [1:size(AineqN,1)];
            lRowsToKeep = ~ismember(lCons,NearOptConstraint(vPrune(1:nObjs)==1));
            AineqN = AineqN(lRowsToKeep,:);
            bineqN = bineqN(lRowsToKeep);
        end
       
        ProbFormNew = ProbForm;
        ProbFormNew.Aineq = AineqN;
        ProbFormNew.Aeq = AeqN;
        ProbFormNew.lb = lbN';
        ProbFormNew.ub = ubN';
        ProbFormNew.bineq = bineqN;
        ProbFormNew.beq = beqN;
    end
 if 0 
    if ~isempty(AMat)
        AMatN = ReorderCols(vKeepDec,AMat);
        
        %Update the right hand side constraints. Substract off from the right-hand side contraint valus the fixed values
        %for the columns to be removed
        [nD lNewDec];
    
        if lNewDec ~= nD
            BrhsN = Brhs - AMat*(vPrune(nObjs+1:n).*vFixedVal(nObjs+1:n)')';
        else
            BrhsN = Brhs;
        end
        
        %Remove rows in the constraints that correspond to near-optimal
        %constraints of removed objectives
        if lNewObjs ~= nObjs
            lCons = [1:size(AMatN,1)];
            lRowsToKeep = ~ismember(lCons,NearOptConstraint(vPrune(1:nObjs)==1));
            AMatN = AMatN(lRowsToKeep,:);
            BrhsN = BrhsN(lRowsToKeep);
        end
        
    else
        AMatN=AMat;
        BrhsN = Brhs;
    end
end
    if ~isempty(cFunc)
        cFuncTemp = ReorderCols(vKeepDec,cFunc);
        %now work on the rows (object functions
        cFuncN = ReorderCols(vKeep(1:lNewObjs),cFuncTemp');
        cFuncN = cFuncN';
    else
        cFuncN = cFunc;
    end
    
    %Work on the BaseAxis
    BaseAxisN = BaseAxis;
    if nObjs>0
        if isempty(find(vKeepObj==BaseAxis(1)))
            BaseAxisN(1)=min(lNewObjs,1);
        else
            BaseAxisN(1) = find(vKeepObj==BaseAxis(1));
        end
    else
        lNewObjs = 0;
    end
    
    if isempty(find(vKeepDec==BaseAxis(2)))
        BaseAxisN(2)=1;
    else
        BaseAxisN(2) = find(vKeepDec==BaseAxis(2));
    end
    
    [mObjs, mDecs]=SplitMatrix(mMatrixN,lNewObjs);
    [vObjLabelsN, vXLabelsN]=SplitMatrix(vLabelsAllN,lNewObjs);
    
    if lNewObjs <= 1
        blShowInset = 0;
    end
        
    %Set parameters that were manipulated to pass via varargin
    varargout = aSetField(varargin_rrax,'tolerance',ToleranceN,'vFixed',vFixedN','vFixedVals',vFixedValN','vStep',vStepN,'vObjLabels',vObjLabelsN','vXLabels',vXLabelsN, ...
                'vXLabelsShort',vXLabelsShortN,'mActCat',mActCatN,'AllowableDeviation',AllowableDeviation,'BaseAxis',BaseAxisN,'mLims',mLimsN,...
                'ProbForm',ProbFormNew,'cFunc',cFuncN,'OptSolRow',OptSolRowN,'NearOptConstraint',NearOptConstraintN,'ShowInsetPlot',blShowInset);
    
    nearoptplotmo2(mObjs, mDecs,varargout{:});
end

function [Indexes] = GetAllValues(hObjs,Prop,Value)
   % find and return the indexes of all the objects in hObjs who property
   % Prop has the value Value. Returns empty if no objects have property
   % with value
   
   %Initialize
   Indexes = [];
   if isempty(hObjs)
       return;
   end
   
   for i=1:length(hObjs)
       try
           cValue = get(hObjs(i),Prop);
           if strcmpi(class(cValue),'numeric')
               if cValue == Value
                   Indexes = [Indexes i];
               end
           else
               if strcmpi(cValue,Value)
                   Indexes = [Indexes i];
               end
           end              
       catch
           continue;
       end
   end
end     

function [Indexes] = GetAllChecked(cbChecks)
    % finds and returns the indexes of all the check boxes in cbChecks
    % array that are checked. If the checkboxes don't exist or the group size is one or all
    % check boxes are unchecked, returns 1 (assumes first box is checked)
    
    %Initialize to first element
    Indexes=1;
    
    if (isempty(cbChecks)) || (length(cbChecks)==1)
        return
    end
        
    for i=2:length(cbChecks)
        if get(cbChecks(i),'Value')==1
            Indexes = [Indexes i];
        end
    end
end

function [Value, Index] = GetFirst(cbChecks,vListValues,IsChecked)
    % finds and returns the first box in the vector of cbChecks
    % boxes that has the value IsChecked (1=checked; 0=not checked). If the checkboxes don't exist or the group size is one, return
    % the first membre of vGroups
    % IsChecked is optional and default value is 1 (checked)
    
        %Find the group to work with
    %OrderedGroups = unique(vGroups);
    %vListValues
    %cbChecks
    if (length(vListValues)==1) || (length(cbChecks)==1) || all((cbChecks==0))
        Value=vListValues(1);
        Index=1;
        return
    end
    
    %Set the value to search for
    if nargin==2 || ~isnumeric(IsChecked)
        cValue = 1; %default to looking for the first checked boxed
    else
        cValue = IsChecked;
    end
    
    %Initialize to first element
    Index=1; %Value = get(cbChecks(Index),'string');
    Value = vListValues(1);
    
    for i=1:length(cbChecks)
        if get(cbChecks(i),'Value')==cValue
            Value = vListValues(1);
            Index=i;
            break
        end
    end
    %if needed, convert to numeric
    if isnumeric(vListValues)
        Value = str2num(Value);
    end
end

function [vGroupSubSet] = ReturnRows(vGroup,vGroupSelect,AllowableDeviationValPU,mMatrix,cbChecks,vFixedVals,sSliders,mConvert,vFixed)
    % Finds and returns a vector of the row numbers in vGroup that have the
    % value vGroupSelect
    %
    % If optional parameters mMatrix, cbChecks, sSliders,
    % are passed, then only returns the rows in mMatrix whose columns have
    % values within AllowableDeviation tolerance for the checked columns.
    %
    % Remember, mMatrix, AllowableDeviation, and vFixedVals are in Data units and sSlider
    % values are in Plot units. Use mConvert to convert from Plot to Data
    % units.
    %

    % OUTPUTS
    % vGroupSubSet = a vector of indexes into vGroup that are the selected
    % rows
    %
    % INPUTS
    % vGruop = an m by 1 vector of groupings
    % vGroupSelect = a parameter or vector of values to search for in vGroup
    % txtAllowableDeviation = handle to a textbox with the error tolerance
    % mMatrix = an m by n matrix of values in data units
    % cbChecks = an n by 1 vector of checkboxes corresponding to the
    % columns in mMatrix. A check mean the column is selected
    % vFixedVals =an n by 1 vector of values in data units
    % vFixed = an 1 by n vector of binaries represented fixed values.
    %   Optional alterantive to cbChecks (use instead)
    
    
    m = size(vGroup,1);
    mInds = [1:m]';
    vGroupSubSetBin = zeros(m,1);
    
    for i=1:length(vGroupSelect)
        if iscell(vGroup)
            vGroupSubSetBin = vGroupSubSetBin+strcmpi(vGroup,vGroupSelect{i});
        else
            vGroupSubSetBin = vGroupSubSetBin+ (vGroup==vGroupSelect(i));
        end
    end
    vGroupSubSet = mInds(vGroupSubSetBin>0);
    
    %AllowableDeviationVal = GetCheckAllowableDeviation(txtAllowableDeviation,'ReturnRows');    
    %AllowableDeviationVal;
    %class(AllowableDeviationVal);
      
    if nargin > 3 %we must also consider the column values within the error tolerance of the fixed values
        %Collect the values from the controls
      
        if nargin < 9
            %use the values from checkboxes
            vFixed = ReadCheckboxValues(cbChecks);
        end
        
        [mF nF] = size(vFixed);
        %vFixedActions = vFixed((nObjs+1):n);
    
      
        [vFixedValsRead vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
        [m2 n2] = size(mMatrix);
    
        if (m2 ~= m) || (n2 ~= length(vFixedVals)) || (n2 ~= length(vFixed))
            warning('nearoptmo2: ReturnRows','Optional parameters have inconsistent dimensions. Ignoring optional inputs')
            return
        end
       
        vFixed;
        vFixedVals;
        
        if sum(vFixed)==0
            return
        end
        
        %Convert AllowableDeviation to Data Units for each column
        AllowableDeviationDU = zeros(1,n2);
        for i=1:n2
            AllowableDeviationDU(i) = ConvertPlot2DataUnits(AllowableDeviationValPU,mConvert(:,i));
        end
        vGroupSubSetWithError = [];
        
        % loop through the rows of mMatrix
        for i=vGroupSubSet'            
               resid = abs(mMatrix(i,:).*vFixed - vFixed.*vFixedVals');
               %[i resid AllowableDeviationValPU]
               
               %resid
               %AllowableDeviationValPU

               if (sum(resid <= AllowableDeviationDU)==n2)
                   vGroupSubSetWithError = [vGroupSubSetWithError; i];
               end
        end
        vGroupSubSet = vGroupSubSetWithError;
    end 
end

function [ProbFormNew,ProbFormExist,AMatNew, BrhsNew, cFuncNew,cFuncFree,vFixed,vFixedVals,AllowableDeviationValue] = UpdateLPMatrix(hWindCurr,nO,lRetType)       
    % Use the checked axes, slider set values, and Subspace Error control setting to define a new linear system of equations representing the feasible area 
    % that is defined by the original linear program components (ProbForm, cFunc) and the current fixed value settings.
    % If an optimal solution and near-optimal tolerance constraint are specified (OptSolRow and NearOptConstraint), 
    % also updates the near-optimal tolerance constraint based on the value in the near-optimal tolerance control and optimal objective function value 
    % 
    % lRetType == A flag to indicate how to update:
    %     0 - subtract fixed values from constraints (reduce matrix)
    %     1 - return fixed values as equity constraints (.Aeq and .beq) in the NewProbForm
    %     2 - Ignore fixed values and leave equity constraints at the original
    %           settings
    %
    % The new system of equations is returned as NewProbForm and
    % cFuncNew
    %
    % They are built as follows:
    % For decision variable axess, this is simply adding two constraints to the AMat:
    %       1) Decision Variable Value <= Fixed Value + Error
    %       2) Decision Variable Value >= Fixed Value - Error
    % For objective functions, also add two constraints to AMat
    %       1) c(X) <= Fixed Value + Error
    %       2) c(X) >= Fixed Value - Error (where c(X) is pulled from the
    %                       coefficiencts in cFunc.
    % For the near-optimal tolerance constraint this is simply
    %       c(X) <= (Optimal Obj Function Value)* (Near Optimal Tolerance)
    %  
    %  As an example:
    %  If the 1st and last decision axes are fixed at FixedValue1 and
    %  FixedValueN, and Objective 1 is fixed at ObjVal1, then the augmented matrix
    %  describing the original and fixed is:
    %           [.Aineq    : .bineq               ]
    %           [1  0 .. 0 : FixedValue1 + Error  ]
    %           [-1 0 .. 0 : -FixedValue1 + Error ]
    %           [0  0 .. 1 : FixedValueN + Error  ]
    %           [0  0.. -1 : -FixedValueN + Error ]
    %           [cFunc     : ObjVal1 + Error      ]
    %           [-cFunc    : -ObjVal1 + Error     ]
    %
    % Note, when Error is zero, the matrix and constraint set is instead
    % reduced and only columns in .Aineq, .Aeq, .lb, and .ub corresponding to unchecked axes are
    % retained.
    %   i.e. [ .Aineq(:,checked==0) : .bineq(:) - .Aineq(:,checked==1)*FixedValue(:)]  
    %    
    % INPUTS
    %  hWindCurr = handle to figure object to pull out stored variables
    %       like AMat, Brhs, cFunc
    %  hControls = cell array of handles to controls on the figure to
    %       pull out current settings like fixed (set) axis, set values,
    %       Error value, etc.
    %  mConvert = matrix of axis conversion factors from plot to data units
    %  nO = number of objectives (columns) in the data set
    %  lReturnType = 1 to include specified equity constraints if
    %       AllowableDeviationValue=0 (otherwise, empty) [default value: 0]
    %
    % OUTPUTS
    %  NewProbForm = updated a-matrix defining the constraint coefficients
    %   (m+2*d x nF) where m=original number of rows/constraints,
    %   d=number of fixed deicision variables + fixed objective function
    %   values, and nF = number of decision variables
    %  BrhsNew = updated right hand side constraint coefficients (m+2*d
    %         x 1)
    %  cFuncNew = nO x nF maxtrix of objective function coefficients for
    %       the nO objectives
    %  cFuncFree = version of cFunc with columns removed for decision
    %       variable values that are fixed (vFixed == 1 and Error=0) when error
    %       value is 0 (reduced mode). Otherwise, this will be the same as
    %       cFuncNew
    %  vFixed = 1 x n vector of binaries with a value of 1 indicating the axis is fixed
    %  vFixedVals = 1 x n vector of values specified the fixed value
    %  AllowableDeviationValue = string value of the Sub Space Error control
    %  AMatEq = d x nF matrix representing decision variables that are
    %       fixed values
    %  BEq = d x 1 vector of fixed variable values

    % Get the handles for the figure controls we will need
    [varargin_fresh,hControls,mObjsCurr] = ReadControls('UpdateLPMatrix',hWindCurr);
    [vFixed, vFixedVals, TickSize, ProbFormExist, cFunc, AllowableDeviationValue, AllowableDeviationValPU, mConvert] = aGetField(varargin_fresh,{'vFixed' 'vFixedVals' 'TickSize' 'ProbForm' 'cFunc' 'AllowableDeviation' 'AllowableDeviationPU' 'mConvert'});   
    [NETolerance, NearOptConstraint,OptSolRow] = aGetField(varargin_fresh,{'Tolerance','NearOptConstraint','OptSolRow'});

    [AMat, Brhs] = OptimiFull(ProbFormExist);
    
    [mF nF] = size(vFixed);
    n = nF;
    %Calc number of decision variables
    nD = nF-nO;
    dInds = [1:nD];
    oInds = [1:nO];

    %Build the new constraint matrix that includes the fixed decision and
    %objective function values
    if isempty(AMat) || isempty(Brhs) || isempty(NearOptConstraint) || isempty(OptSolRow)
        drawnow;
        msgbox('Need additional input data to generate new alternatives with Matlab. Specify ProbForm (to define inequity, equity, lower-bound, and/or upper-bound constraints), NearOptConstraint (row in inequity constraints that is the near-optimal constraint), and OptSolRow (row in original mDecisions that is the optimal solution) on opening the plotter. Instead, set Generate engine to Data or GAMS.','Error');
        AMatNew = [];
        ProbFormNew = [];
        return
    end

    AMatNew = ProbFormExist.Aineq;
    BrhsNew = ProbFormExist.bineq;
    cFuncNew = cFunc;
    cFuncFree = cFunc;
    
    AMatEq = ProbFormExist.Aeq;
    BEq = ProbFormExist.beq;
    lB = ProbFormExist.lb;
    uB = ProbFormExist.ub;
    
    %Determine lMatRet
    if nargin < 3
        lRetType = 0;
    end

    %Update the near-optimal tolerance constraint(s)
    if (nO > 0) && ~isempty(NearOptConstraint) && ~isempty(OptSolRow)
        %Error checking on Near-optimal tolerance
        Tolerance = abs(NETolerance); %Screen out negative values
        for i=1:nO
            %Loop through objectives
            if Tolerance(i) < 1
                %We have a maximization function; add a negative sign to
                %preserve the integrety of the <= constraint!
                Tolerance(i) = -Tolerance(i);
            end
        
            BrhsNew(NearOptConstraint(i)) = Tolerance(i)*mObjsCurr(OptSolRow(i),i);
        end
    end   
    
    dToAdd = dInds(vFixed(nO+1:n) ~= 0);
    oToAdd = oInds(vFixed(1:nO) ~= 0);
    
    if AllowableDeviationValPU==0
       % fixed values for decision variables
       if lRetType == 0
           %Reduce matrix -- we need to remove the variable from the constraint set and use a
            %modified (reduced) matrix set
           vFixedValsDs = vFixedVals(nO+1:end);
           BrhsNew = BrhsNew - AMatNew(:,vFixed(nO+1:end)==1)*vFixedValsDs(vFixed(nO+1:end)==1);
           AMatNew = AMatNew(:,vFixed(nO+1:end)==0);
           if ~isempty(lB)
               lB = lB(vFixed(nO+1:end)==0);
           end
           if ~isempty(uB)
               uB = uB(vFixed(nO+1:end)==0);
           end           
           if ~isempty(cFunc)
               cFuncFree = cFunc(vFixed(1:nO)==0,vFixed(nO+1:end)==0);
           end
       elseif lRetType == 1
          %Create equity constraints
          BEq = vFixedVals(nO+dToAdd);
          for i=1:length(dToAdd)
              AMatEq = [AMatEq; circshift(eye(nD,1),dToAdd(i)-1)'];
          end
       end
    else
        %Add constraints to allow decision variables within some
        %error/tolerance of the fixed values
        blIsReduced = 0;
        for i=dToAdd
            cRow = circshift(eye(nD,1),i-1)';
            ErrorValDU = ConvertPlot2DataUnits(AllowableDeviationValPU,mConvert(:,nO+i));
            vFixedVals(i+nO);
            AMatNew = [AMatNew;cRow;-cRow];
            BrhsNew = [BrhsNew; vFixedVals(i+nO) + ErrorValDU; - vFixedVals(i+nO) + ErrorValDU];
        end
    end

    if (sum(vFixed(1:nO)) > 0) && (isempty(cFunc) || size(cFunc,1) < sum(vFixed(1:nO)))
        warning('Need to specify parameter cFunc to fix objective function values or not enough cFuncs specified. Continuing without')
    else
        for i=oToAdd
            ErrorValDU = ConvertPlot2DataUnits(AllowableDeviationValPU,mConvert(:,i));
            AMatNew = [AMatNew;cFuncNew(i,:);-cFuncNew(i,:)];
            BrhsNew = [BrhsNew; vFixedVals(i) + ErrorValDU; - vFixedVals(i) + ErrorValDU];     
        end
    end
    
   %Move results back into the new ProbForm structure
   ProbFormNew = ProbFormExist;
   ProbFormNew.Aineq = AMatNew;
   ProbFormNew.bineq = BrhsNew;
   ProbFormNew.Aeq = AMatEq;
   ProbFormNew.beq = BEq;
   ProbFormNew.lb = lB;
   ProbFormNew.ub = uB;    
end

function [mResCompact,mResultsData,cRowsRet] = ReturnDataExtents(hCurr,mData,nO,blStrictExtent)
     %Query the max and min values from the data records in mData
     % defined by cbChecks (fixed axes) within AllowableDeviation of
     %vFixedVals using ReturnRows. Find a min and max for each axis. returns 2*n records
     %
     % blStrictExtent = 1 will return the strict extents along each axis,
     %       i.e., values for fixed axes at their fixed values
     % blStrictExtent = 0 will unfix each fixed axis in turn and return the
     %      extents for that axis (with other fixed axes still fixed)

     % mResCompact is the results in plot units in compact form, i.e., 2 x n matrix
     %     with min in row 1 and max in row 2
     % mResultsData is full results in a 2n x n matrix 
     %     with rows and 1 and 2 the data rows containing the mins and maxes for axis 1, rows 3
     %     and 4 the data rows containing the min and max for variable 2, etc.
     % cRows = row indexes into the original data set that are returned
     % blStrictExtent = a binary variable with a value of 1 indicates
     %     strict extents and 0 indicates fixed variables are related in turn

     % Read in the values from controls
   
    [varargin_fresh,hControls] = ReadControls('ReturnDataExtents',hCurr);    
    [AllowableDeviationVal,AllowableDeviationValPU,mGroupDataSort,vGroup,vFixedVals,mConvert] = aGetField(varargin_fresh,{'AllowableDeviation','AllowableDeviationPU','mGroupData','vGroup','vFixedVals','mConvert'});
    
    %Read in select controls
    [cbChecks,sSlider] = aGetField(hControls,{'AxisChecked','Sliders'});
    
    [m,n] = size(mData);
    
    %[AllowableDeviationVal, AllowableDeviationValPU] = GetCheckAllowableDeviation(txtAllowableDeviation,'ReturnDataExtents',TickSize);
    ColInds = [1:n];

    %Find the groups to work with
    selGroups = mGroupDataSort(cell2mat(mGroupDataSort(:,2))==1,1); %unique(vGroup);
    cRows = ReturnRows(vGroup,selGroups,AllowableDeviationValPU,mData,cbChecks,vFixedVals,sSlider,mConvert);
    %mResultsData = mData(cRows,:);
    %size(mResults)

    %Build the maximum extents for each axis. For unfixed variables
    %these extents are the same extents previously identified. For
    %fixed variables, we unfix that particular variable and resolve

    %Cull all for the unfixed variables
    if length(cRows)==1
        mResCompact = [mData(cRows,:);mData(cRows,:)];
        cRowsRet = cRows;
    else
        [mDataMin,cIMin] = min(mData(cRows,:));
        [mDataMax,cIMax] = max(mData(cRows,:));
        mResCompact = [mDataMin;mDataMax];
        cRowsRet = [cRows(cIMin)';cRows(cIMax)'];
    end
        
    if blStrictExtent==0
        %loop through the fixed variables and unfix each in turn to obtain
        %wider extent
        for i=ColInds(vFixed==1)
            vFixedUse = vFixed;
            vFixedUse(i) = 0;
            cRowsCurr = ReturnRows(vGroup,selGroups,AllowableDeviationValPU,mData,cbChecks,vFixedVals,sSlider,mConvert,vFixedUse);
            if length(cRowsCurr)==1
                mResCompact(:,i) = repmat(mData(cRowsCurr,i),2,1);
                cRowsRet(:,i) = repmat(cRowsCurr,2,1);
            else
                [cMin, iMin] = min(mData(cRowsCurr,i));
                [cMax, iMax] = max(mData(cRowsCurr,i));
                mResCompact(:,i) = [cMin;cMax];
                cRowsRet(:,i) = [iMin;iMax];
            end
        end
    end
    
    %Reshape cRowsRet into a row vector
    if (~isempty(cRowsRet)) && (size(cRowsRet,1) > 1)
        cRowsRet = reshape(cRowsRet,2*n,1);
    end
    
    mResultsData = mData(cRowsRet,:);
end

%function Resample(hWind,event,txtTolerance,txtAllowableDeviation,txtNumSamples,cbGroupChecks,cbChecks,sSliders,mConvert,mData,sControls,nO,hWindCurr)
function [mObjsNew, mDecsNew] = Resample(hWind,event,hWindCurr,mData,nO)
    %Re-sample an additional specified number of alternatives from the sub-space defined by the checked axes and fixed values associated with those
    %axes. For decision variable axess, this is simply adding two constraints to the AMat:
    %       1) Decision Variable Value <= Fixed Value + Error
    %       2) Decision Variable Value >= Fixed Value - Error
    % For objective functions, also add two constraints to AMat
    %       1) c(X) <= Fixed Value + Error
    %       2) c(X) >= Fixed Value - Error (where c(X) is pulled from the
    %                       coefficiencts in cFunc.
    % Then run the random sample routine.
    %  If the 1st and last decision axes are fixed at FixedValue1 and
    %  FixedValueN, and Objective 1 is fixed at ObjVal1, then the augmented matrix
    %  describing the original and fixed is:
    %           [AMat      : Brhs                 ]
    %           [1  0 .. 0 : FixedValue1 + Error  ]
    %           [-1 0 .. 0 : -FixedValue1 + Error ]
    %           [0  0 .. 1 : FixedValueN + Error  ]
    %           [0  0.. -1 : -FixedValueN + Error ]
    %           [cFunc     : ObjVal1 + Error      ]
    %           [-cFunc    : -ObjVal1 + Error     ]
    %
    % Note, when Error is zero, the matrix and constraint set is instead
    % reduced and only columns in AMat corresponding to unchecked axes are
    % retained.
    %   i.e. [ AMat(:,checked==0) : Brhs(:) - AMat(:,checked==1)*FixedValue(:)]  
    %    
    % INPUTS
    %  hWindCurr = handle to figure object to pull out stored variables
    %  hControls = cell array of handles to controls on the figure
    %  mConvert = matrix of axis conversion factors from plot to data units
    %  mData = the data in data units
    %  nO = number of objectives (columns) in the data set
                        
    % Get the handles for the figure controls we will need
    
    %Reading in the controls and input data
    [varargin_fresh,hControls] = ReadControls('Resample',hWindCurr);
    %Update the constraint matrix
    [ProbNew,ProbOrig,AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,AllowableDeviationValue] = UpdateLPMatrix(hWindCurr,nO);

    if isempty(AMatNew)
        %There was an error
        return
    end
    
    [ToleranceValue, sGamsFile, NumSamples, ObjSamplesPercent, vGroup, mActCat, mGroupDataSort,AllowableDeviationValPU, ErrorResid] = aGetField(varargin_fresh,{'Tolerance' 'sGamsFile' 'NumSamples' 'ObjSamplesPercent' 'vGroup' 'mActCat' 'mGroupData' 'AllowableDeviationPU' 'ErrorResid'});
    
    [m n] = size(mData);
    nD = n-nO;
    
    %Split the samples among objective (specified percentage) and decision
    %variable axes
    
    ObjDecSplit = NumSamples*[100-ObjSamplesPercent ObjSamplesPercent]/100; 
            
    %Check extents of matrixes
    [mExtOrig, mExtCompactOrig] = maxextentind(ProbOrig);
    [mExtNew,mExtCompactNew] = maxextentind(ProbNew);
      
    %Resample
    [NewSols,vValid] = stratgibbs(ObjDecSplit,ProbNew,struct('matformat','reduce','lincombo',cFuncFree','MCMCMethod','cprnd-gibbs','errorresid',ErrorResid));
    %[NewSols] = cprnd(NumSamples,AMatNew,BrhsNew,struct('method','hitandrun'));
    %vValid = ones(NumSamples,1);

    %Only use the valid ones, screen out rows
    NewSols = NewSols(vValid==1,:);
    mDecsNew = NewSols;
    [mNS nNS] = size(NewSols);
    
    if mNS==0
        mObjsNew = [];
        mDecsNew = [];
        return;
    end
    
    %B. Work on the decision variables
    if (nNS < nD) && (nNS + sum(vFixed(nO+1:end)) == nD)
        mExtComNF = CompactToFull(mExtCompactNew',vFixed(nO+1:end),vFixedVals(nO+1:end)');
        mExtCompactNew = mExtComNF';
        mDecsNew = CompactToFull(NewSols,vFixed(nO+1:end),vFixedVals(nO+1:end)');    
    end
       
    %Print out new decision ranges
    %sprintf('New Decision Ranges')
    ['# Location Orig Extents  New Extents']
    %[[1:nD]' mExtCompactOrig mExtComNF]
    [num2cell([1:nD]') mActCat num2cell(mExtCompactOrig) num2cell(mExtCompactNew)]
    
    %C. Work on the Objective functions; have to calculate objective function values for the new solutions
    if nO>0 && ~isempty(cFunc)
        %[size(mDecsNew); size(cFunc')]
        mObjsNew = mDecsNew*cFunc';       
        %Print out the range of new sampled obj function values
        sprintf('New Objective Function Range')
        [min(mObjsNew); max(mObjsNew)]
        %figA = figure;
        %hist(mObjsNew(m+1:end))
        %hold on
        %title('Histogram of sampled objective function values')
        %hold off
    else
        mObjsNew = [];
    end 
end 

function [mObjsNew, mDecsNew] = EnumerateWithGams(hWind,event,hWindCurr,mData,nO)
    %Generate additional solutions with the enumeration algorithm that
    %calls the specified GAMS file
    %    
    % INPUTS
    %  hWindCurr = handle to figure object to pull out stored variables
    %  hControls = cell array of handles to controls on the figure
    %  mConvert = matrix of axis conversion factors from plot to data units
    %  mData = the data in data units
    %  nO = number of objectives (columns) in the data set
    
    % OUTPUTS
    % mObjNew = m x nO matrix of new objective function values generated by
    %    GAMS in data units. Each row m is a new solution
    % mDecNew = m x n matrix of decision variable values generated by GAMS
    %    in data units. Each row m is a new solution.
    %
                        
    % Get the handles for the figure controls we will need
    [varargin_fresh,hControls] = ReadControls('EnumerateWithGams',hWindCurr);
    [ToleranceValue, sGamsFile, NumSamples, vFixed, mConvert, lRunMode] = aGetField(varargin_fresh,{'tolerance' 'sGamsFile' 'NumSamples' 'vFixed' 'mConvert' 'GenerateType'});
    %[txtTolerance txtNumSamples, cbChecks, txtAllowableDeviation, sSliders, txtGroupOrders,txtGroupNames,cbGroupChecks,txtGamsFile] = aGetField(hControls,{'Tolerance','NumSamples','AxisChecked','AllowableDeviation','Sliders','GroupOrders','GroupNames','GroupChecks','GamsFile'});
    
    [mF nF] = size(vFixed);
    [m n] = size(mData);
    
    %Calc number of decision variables
    nD = n-nO;
    dInds = [1:nD];
    oInds = [1:nO];
        
    %Read in the variable arguments
    %varargin_fresh = getappdata(hWindCurr,'varargs');
    [vFixedVals,vGroup,TickSize,vXLabelsShort,mGroupDataSort,AllowableDeviationValue] = aGetField(varargin_fresh,{'vFixedVals' 'vGroup' 'TickSize' 'vXLabelsShort' 'mGroupData' 'AllowableDeviation'});   
    
    %Update the groups
    %mGroupDataSort = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
    
    %[AllowableDeviationValue, AllowableDeviationValPU] = GetCheckAllowableDeviation(txtAllowableDeviation,'Resample',TickSize);
    %[vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    %Change the directory/folder to the one the GAMS File is in
    if exist(sGamsFile,'file') == 0
        msgbox(['The file ',sGamsFile,' does not exist. Cancelling.'],'Warning');
        return
    end
    [sPath,sFileName,sFileExt] = fileparts(sGamsFile);
    if ~strcmpi(sPath,'')
        cd(sPath);
    end
    
    %[vXLabelsShort' num2cell(vFixed(nO+1:end)') num2cell(vFixedVals(nO+1:end))]
 
    %Pass the values and enumerate the integer solutions
    [mObjs, mResultsInt, mResultsVal, uelsOut, ReturnFlag, GamsStats, NumSolvs] = EnumNEIntSolsGams4(sGamsFile,vFixed(nO+1:end),vFixedVals(nO+1:end),vXLabelsShort,ToleranceValue,2,lRunMode);
    
    %[min(mResultsVal);max(mResultsVal)]
    
    %[vObjs mResultsVal ReturnFlag]
    
    [mNSolv nNSolv] = size(mResultsVal); %Number of solves
           
    mDecsNew =mResultsVal(ReturnFlag>0,:);  
    mObjsUse = mObjs(ReturnFlag>0);

    %C. Work on the Objective functions; have to calculate objective function values for the new solutions
    if nO>0 && ~isempty(mObjsUse)
        mObjsNew = [1.41*mObjsUse];
    else
        mObjsNew = mObjs;
    end
       
    %Show errors, if any
    numErrors = sum(ReturnFlag<0);

    %[mObjsNew mDecisionsNew]
    [sprintf('%d solves made\n%d new alternatives added\n%d errors',NumSolvs,NumSolvs-numErrors,numErrors)]

    if numErrors>0
        ['Alternatives with errors']
        [ReturnFlag(ReturnFlag<0) GamsStats(ReturnFlag<0,:) mResultsInt(ReturnFlag<0,:)]
    end        
end  


function GenerateNewAlts(hWind,event,hWindCurr)
    % Generates new alternatives based on the selctions in the Generate
    % New Alternatives Box on the Interact Tab. Adds as a new group to the plot.

    %Read the controls
    [varargin_fresh,hControls,mObjs,mDecisions] = ReadControls('GenerateNewSols',hWindCurr);    
    [vGroup,lGenUse,lGenType,mConvert,nO] = aGetField(varargin_fresh,{'vGroup' 'GenerateMethod' 'GenerateType' 'mConvert' 'nO'});   
    mData = [mObjs mDecisions];
    
    cTypeSuffix = {'Sing' 'Ext' 'Resample' 'Enum'};
    cUseSuffix = {'Dat' 'Mat' 'Gams'};
    
    cTypeSuffixFull = {'Single Solution' 'Maximum Extents' 'Resample' 'Eunumerate'};
    cUseSuffixFull = {'Data' 'Matlab' 'Gams'};
    
    if lGenUse==1
        %Data
        if lGenType~=2
            msgbox(['No function defined for GenType #', num2str(lGenType), '. Ignoring'],'Warning');
            return;
        else
           [mResCompact,mResultsData, cRows] = ReturnDataExtents(hWindCurr,mData,nO,1);
           [mObjNew,mDecNew] = SplitMatrix(mResultsData,nO);
        end
    elseif lGenUse==2
        %Matlab functions calling on the definition of the linear model in
        %the parameters AMat, Brhs, and cFunc.
        
        %First update the definition of the contraints based on the control
        %settings (checkboxes, fixed values, etc.
        [ProbNew,ProbOrig,AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,AllowableDeviationValue] = UpdateLPMatrix(hWindCurr,nO);
        
        if isempty(AMatNew)
            %There was an error
            return
        end
       
        if lGenType==1
            %Single solution, optimize
            %[mDecNew, mObjNew, exitflag, minoutput] = linprog(cFunc,AMatNew,BrhsNew,[],[],[],[],[],struct('maxiter',1000,'Display', 'off'));
            ProbLPNew = ProbNew;
            ProbLPNew.f = cFunc;
            [mDecNew, mObjNew, exitflag, minoutput] = linprog(ProbLPNew);
            
        elseif lGenType==2
            %Maximum extents   
            [mDecNew,mResCompact] = maxextentind(ProbNew);
            
            %Calculate objective function values for new solutions
            [mObjNew] = mDecNew*cFunc';

        elseif lGenType==3
            %Resample
            [mObjNew, mDecNew] = Resample(0,0,hWindCurr,mData,nO);
        else
            %Enumerate
            warning(['No function defined for GenType #', num2str(lGenType), '. Ignoring']);
            return;   
        end
    elseif lGenUse==3
        %GAMS
        if lGenType==3
            warning(['No function defined for GenType #', num2str(lGenType), '. Ignoring']);
            return;
        else
            [mObjNew, mDecNew] = EnumerateWithGams(0,0,hWindCurr,mData,nO);

        end
    end
        
    sSuffix = [cTypeSuffix{lGenType},cUseSuffix{lGenUse}];
    
    %Combine the new solutions with existing ones and add to the plot as a new group
    
    [mNS nNS] = size(mDecNew);
    
    if mNS >= 1
        %We have some new solutions
        
        if ismember(lGenType, [1 2])
            %Show the new solutions
            [mObjNew mDecNew]
        end

        %A. Work on the decision variables
        %Map the output columns to the input columns
        %vMap = MapLabels(uelsOut,vXLabelsShort);

        mDecisionsAll = [mDecisions; mDecNew];
        mObjsAll = [mObjs; mObjNew];
        %Print out results in the order the decision variables appear on
        %the axes
        %[ReturnFlag GamsStats mResultsInt]
        %[ReturnFlag mResultsVal(:,vMap)]

        %Print out new decision ranges
        %sprintf('New Decision Ranges')
                %D. Work on the grouping. The new group will appear with a child name just below the
        %first checked group in the resorted list

        %Find the first checked group
        [sNewGroupName, mGroupDataNewSort] = GenerateNewGroup(hWindCurr,'',sSuffix);

        %Assign the new group name to the enumerated records
        if iscell(vGroup)
            vGroupAdd = repmat({sNewGroupName},mNS,1);        
        else
            vGroupAdd = sNewGroupName*ones(mNS,1);
        end
        
        sprintf('%s using %s: %d solutions added',cTypeSuffixFull{lGenType},cUseSuffixFull{lGenUse},mNS)
        % Update the variables
        varargout = aSetField(varargin_fresh,'vGroup',[vGroup;vGroupAdd],'mGroupData',mGroupDataNewSort(:,1:3));  

        nearoptplotmo2(mObjsAll,mDecisionsAll,varargout{:});
    else
        h = errordlg('No new alternatives generated','Solution Error');

        uiwait(h); 
    end
end

function AddObjectives(hind,evnt,hCurr,dfAns)
    % Solicits the user for new objective function(s), updates the model
    % formulation, and resamples alternatives according to the updated
    % formulation. There must be an existing model formulation defined by
    % AMat and Brhs.
    
    % INPUTS
    %   hCurr = handle to the current figure window
    %   dfAns = optional cell array of strings with default answers
    %   (default is empty)
    
    [VarsRet, hCont, mObjsOrig, mDecsOrig] = ReadControls('AddObjectives',hCurr);
    [ProbForm,cFuncOrig,NearOptConstraintOrig,OptSolRowOrig,vObjLabelsOrig, vToleranceOrig, sAxisScales, mLimsOrig, nOOrig, vGroups, blShowInset] = aGetField(VarsRet, {'ProbForm' 'cFunc'  'NearOptConstraint' 'OptSolRow' 'vObjLabels' 'tolerance', 'AxisScales' 'mLims' 'nO' 'vGroup' 'ShowInsetPlot'});
    
    %[lCon,nD] = size(AMatOrig);
    if isempty(ProbForm)
        warning('Need problem data for exsiting problem (ProbForm specifying constraints) to add objectives')
    else
        [lCon,nD] = size(ProbForm.Aineq);
    end
    
    %Step 1. Solicit the user for updated model information
    lNumFields = 6;
    options.interpretor = 'tex';    
       
    %Fields are:
    %   1 - number of new objectives (nNew)
    %   2 - Objective function coefficients (nNew x nD)
    %   3 - Labels (1 x nNew)
    %   4 - Units of measurement (1 x nNew)
    %   5 - Direction (1 x nNew)
    %   6 - Tolerance (1 x nNew)
    %   7 - Plot limits (2 x nNew) (only if AxisScales is custom)
    
    cPrompts = {['Enter information for one or more new objectives (nNew). All entries must be valid Matlab commands and can reference data, functions, or existing variables/matrices in the Base Workspace (note: ',num2str(nD),' = number of decision variables)', ...
        char(10), char(10) 'Number of new objectives (nNew; integer):'], ...
        ['New objective function coefficients (nNew x ',num2str(nD),'):'], ...
        'Labels (1 x nNew cell array of strings; e.g., {''label 1'' ''label 2''}):', ...
        'Units of measurement (1 x nNew cell array of strings):', ...
        ['Direction (min or max) (1 x nNew cell array of string):'], ...
        ['Near-optimal toleraence (fraction) (1 x nNew):']};
    
    if strcmpi(sAxisScales,'custom')
        %Add an additional input for the custom limits
        cPrompts = {cPrompts{:}, 'Plot limits (2 x nNew):'};
        lNumFields = lNumFields+1;
    end
   
    if (nargin <= 3) || (length(dfAns) ~= lNumFields)
        dfAns = repmat({''},lNumFields,1);
    end
    drawnow; pause(0.05); % Show a clean dialog
    vAnswers = inputdlg(cPrompts,'Add New Objective Functions',1,dfAns,options);
    
    %Step 2. Check validity of input
     if isempty(vAnswers) 
        %Canceled, exit immediately
        return
     elseif all(strcmpi(vAnswers,''))
        %No input provided, warn
        h = errordlg('No imput provided on new objective functions. Canceling','Input Error');
        uiwait(h);
        return
    end
    
    %Go though each input and check validity
    cAnswers = cell(lNumFields,1);
    mSizes = zeros(lNumFields,2);
    for i=1:lNumFields
       try 
          cAnswers{i} = eval(vAnswers{i}); 
       catch
          h = errordlg(['Error with input: ',vAnswers{i},'. Try again.'],'Input Error');
          uiwait(h);
          vAnswers{i} = '';
          AddObjectives(0,0,hCurr,vAnswers);
          return
       end
       mSizes(i,:) = size(cAnswers{i});
    end
    
    %Check for consistency with existing model inputs
    nNewO = cAnswers{1}; %number of new objectives
    lFields = [2 2:lNumFields];
    lTrueSize = [nD repmat(nNewO,1,lNumFields-2+1)];
    iSzCol = [2 1 2*ones(1,lNumFields-2)];
    for i=lFields
        if (lTrueSize(i) ~= mSizes(i,iSzCol(i)))
            sDir = 'rows';
            if iSzCol(i) == 2
                sDir = 'columns';
            end
            
            h = errordlg(['Error with input: ',vAnswers{i},'. Number of ', sDir, ' must be ',num2str(lTrueSize(i)),'. Try again.'],'Input Error');
            uiwait(h);
        	vAnswers{i} = '';
            AddObjectives(0,0,hCurr,vAnswers);
            return
        end
    end      
    
    %Check consistency of direction and tolerance and calculate optimal
    %objective function value
    vOptVal = zeros(1,nNewO);
    dirs = zeros(1,nNewO);
    vOptSols = [];
    
    ProbForm.solver = 'linprog';
    ProbForm.options = struct('Algorithm','interior-point','maxiter',15000,'Display','off');
    dirs = 2*strcmpi(cAnswers{5},'min')-1; % Map boolean (1 0) to (1 -1)
 
    for i=1:nNewO
             
        if (dirs(i)*cAnswers{6}(i) < dirs(i)*1)
            warning(['Tolerance of new objective #', num2str(i),' does not match direction. Inverting tolerance'])
            cAnswers{6}(i) = 1/cAnswers{6}(i);
        end
        
        ProbFormSingle = ProbForm;
        ProbFormSingle.f = dirs(i)*cAnswers{2}(i,:);
        [a,b,c,d] = linprog(ProbFormSingle);
        %[a,b,c,d] = linprog(dirs(i)*cAnswers{2}(i,:),AMatOrig,BrhsOrig);
        vOptSols = [vOptSols; a'];
    end
    
    %Step 3. Update the model formulation    
    %Objective function labels
    vObjLabelsAdd = cell(1,cAnswers{1});
    
    for i=1:nNewO
        vObjLabelsAdd{i} = strcat(cAnswers{3}(i),' (',cAnswers{4}(i),')');
    end
    
    vObjLabelsNew = [vObjLabelsOrig vObjLabelsAdd{:}];
    
    %Objective function coefficients, Near-optimal tolerable deviation
    %constraints, tolerance values, and axis limits
    mNewObjCoef = cAnswers{2};
    vOptObjVals = mNewObjCoef*vOptSols';

    cFuncNew = [cFuncOrig; mNewObjCoef];
    ProbFormNew = ProbForm;
    ProbFormNew.Aineq = [ProbFormNew.Aineq; mNewObjCoef.*repmat(dirs',1,nD)];
    ProbFormNew.bineq = [ProbFormNew.bineq; diag(vOptObjVals).*dirs'];
    NearOptConstraintNew = [NearOptConstraintOrig length(ProbForm.bineq)+[1:nNewO]];
    %AMatNew = [AMatOrig; mNewObjCoef.*repmat(dirs',1,nD)];
    %BrhsNew = [BrhsOrig; diag(vOptObjVals).*dirs'];
    %NearOptConstraintNew = [NearOptConstraintOrig length(BrhsOrig)+[1:nNewO]];
    OptSolRowNew = [OptSolRowOrig size(mObjsOrig,1)+[1:nNewO]];
    ToleranceNew = [vToleranceOrig cAnswers{6}];
    
    %Assign a new group name to the new solutions and update the groupings
    [sNewGroupName, mGroupDataNewSort] = GenerateNewGroup(hCurr,'NewObj');
    vGroupNew = [vGroups; repmat({sNewGroupName},nNewO,1)];

    if strcmpi(sAxisScales,'custom')
        [mLimsObj,mLimsDec] = SplitMatrix(mLimsOrig,nOOrig);
        mLimsNew = [mLimsObj cAnswers{7} mLimsDec];
    else
        mLimsNew = mLimsOrig;
    end
    
    %Insert new default values for new objectives
    [vFixedOrig,vFixedValsOrig,vStepOrig] = aGetField(VarsRet, {'vFixed' 'vFixedVals' 'vStep'});
    vParamsToSplit = {vFixedOrig,vFixedValsOrig',vStepOrig};
    vParmsToInsert = {0 0 1/20};
    vParamsNew = cell(length(vParamsToSplit),1);
    
    for i=1:length(vParamsToSplit)
        [vObjs,vDecs] = SplitMatrix(vParamsToSplit{i},nOOrig);
        vParamsNew{i} = [vObjs vParmsToInsert{i} vDecs];
    end    
    
    %Update the data set
    mObjsNew = [mObjsOrig mDecsOrig*mNewObjCoef'; vOptSols*cFuncOrig' vOptObjVals]; %calculate new objective function values for prior alternatives; include objective function values for new optimal solutions
    mDecsNew = [mDecsOrig; vOptSols];
    
    %Save the new parameters and dataset back to the Window so we can use them in resampling
    VarsNew = aSetField(VarsRet,'ProbForm',ProbFormNew,'cFunc',cFuncNew,'NearOptConstraint',NearOptConstraintNew,'OptSolRow',OptSolRowNew, ... %'AMat',AMatNew,'Brhs',BrhsNew
            'vObjLabels', vObjLabelsNew,'Tolerance',ToleranceNew,'mLims',mLimsNew,'nO',nOOrig+nNewO,'vGroup',vGroupNew,'mGroupData',mGroupDataNewSort(:,1:3), ...
            'vFixed', vParamsNew{1},'vFixedVals',vParamsNew{2}','vStep',vParamsNew{3});       
    setappdata(hCurr,'varargs',VarsNew);
    setappdata(hCurr,'mObjs',mObjsNew);
    setappdata(hCurr,'mDecs',mDecsNew);   
    
    %Step 4. Rsample and Replot
    GenerateNewAlts(0,0,hCurr);
    %ToggleInset(aGetField(hCont,{'ShowInsetPlot'}),0); %Show the inset plot
end    

function AddConstraints(hind,evnt,hCurr,dfAns)
    % Solicits the user for new constraint(ss) in the form of A x <=,==,>= b, updates the model
    % formulation, and resamples alternatives according to the updated
    % formulation. There must be an existing model formulation defined by
    % ProbForm.
    
    % INPUTS
    %   hCurr = handle to the current figure window
    %   dfAns = optional cell array of strings with default answers
    %   (default is empty)
    
    [VarsRet, hCont] = ReadControls('AddConstraints',hCurr);
    [ProbForm] = aGetField(VarsRet, {'ProbForm'});
    %[AMatOrig,BrhsOrig] = aGetField(VarsRet, {'AMat' 'Brhs' });
    
    [AMat,Brhs] = OptimiFull(ProbForm);
    
    [nD] = size(AMat,2);
    
    %Step 1. Solicit the user for updated model information
    lNumFields = 5;
    options.interpretor = 'tex';    
       
    %Fields are:
    %   1 - number of new contraints (nNew)
    %   2 - Left hand side coefficiencts (A)(nNew x nD)
    %   3 - Right hand side coefficients (B) (nNew x 1)
    %   4 - Direction (nNew x 1)
    %   5 - Group name
    
    cPrompts = {['Enter information for one or more new constraints (nNew). Constraints are in the form Ax <=> b. All entries must be valid Matlab commands and can reference data, functions, or existing variables/matrices in the Base Workspace (note: ',num2str(nD),' = number of decision variables)', ...
        char(10), char(10) 'Number of new constraints (nNew; integer):'], ...
        ['New left-hand side coefficients (A) (nNew x ',num2str(nD),'):'], ...
        ['New right-hand side coefficients (b) (nNew x 1):'], ...
        ['Direction (-1, 0, or 1 to indicate >=, =, or <=) (1 x nNew):'], ...
        ['Group name:']};
    
   
    if (nargin <= 3) || (length(dfAns) ~= lNumFields)
        dfAns = repmat({''},lNumFields,1);
    end
    drawnow; pause(0.05); % Show a clean dialog
    vAnswers = inputdlg(cPrompts,'Add New Constraints',1,dfAns,options);
    
    %Step 2. Check validity of input
     if isempty(vAnswers) 
        %Canceled, exit immediately
        return
     elseif all(strcmpi(vAnswers,''))
        %No input provided, warn
        h = errordlg('No imput provided on new constraints. Canceling','Input Error');
        uiwait(h);
        return
    end
    
    %Go though each input and check validity
    cAnswers = cell(lNumFields,1);
    mSizes = zeros(lNumFields,2);
    for i=1:lNumFields
       try 
          cAnswers{i} = eval(vAnswers{i}); 
       catch
          h = errordlg(['Error with input: ',vAnswers{i},'. Try again.'],'Input Error');
          uiwait(h);
          vAnswers{i} = '';
          AddConstraints(0,0,hCurr,vAnswers);
          return
       end
       mSizes(i,:) = size(cAnswers{i});
    end
    
    %Check for consistency with existing model inputs
    nNewC = cAnswers{1}; %number of new objectives
    lFields = [2 2:lNumFields-1];
    lTrueSize = [nD repmat(nNewC,1,lNumFields-2)];
    iSzCol = [2 ones(1,lNumFields-2)];
    for i=lFields
        if (lTrueSize(i) ~= mSizes(i,iSzCol(i)))
            sDir = 'rows';
            if iSzCol(i) == 2
                sDir = 'columns';
            end
            
            h = errordlg(['Error with input: ',vAnswers{i},'. Number of ', sDir, ' must be ',num2str(lTrueSize(i)),'. Try again.'],'Input Error');
            uiwait(h);
        	vAnswers{i} = '';
            AddConstraints(0,0,hCurr,vAnswers);
            return
        end
    end      
        
    %Step 3. Update the model formulation
    ProbFormNew = ProbForm;
   
    %Loop through each contrainst
    for l=1:cAnswers{1}
        if cAnswers{4}(l) == 0
            % equity constraint, update the equity constraint parameters
            ProbFormNew.Aeq = [ProbFormNew.Aeq; cAnswers{2}(l,:)];
            ProbFormNew.beq = [ProbFormNew.beq; cAnswers{3}(l)];            
        else
            %less than or equal to or greater than or equal to constraint
            %Direction indicated by Direction variable
            ProbFormNew.Aineq = [ProbFormNew.Aineq; cAnswers{4}(l)*cAnswers{2}(l,:)];
            ProbFormNew.bineq = [ProbFormNew.bineq; cAnswers{4}(l)*cAnswers{3}(l)];     
        end
    end
      
    %Save the new parameters and dataset back to the Window so we can use them in resampling
    VarsNew = aSetField(VarsRet,'ProbForm',ProbFormNew); %'AMat',AMatNew,'Brhs',BrhsNew);       
    setappdata(hCurr,'varargs',VarsNew);
    
    %Step 4. Rsample and Replot
    GenerateNewAlts(0,0,hCurr);
end 

function SaveFig(hind,evnt,hCurr)
    % Save the data (mObjs, mDecisions, varargs) in the figure with handle hCurr to the base workspace
    % Reads controls and generates the current variable arguement in
    % (VarsRet). Similarly retreives the mObjectives and mDecisions
    %
    % INPUTS
    % hCurr = handle of specified figure
    % 
    
    
    [VarsRet, hCont] = ReadControls('SaveFig',hCurr);
    setappdata(hCurr,'varargs',VarsRet);
    mObjsRet = getappdata(hCurr,'mObjs');
    mDecRet = getappdata(hCurr,'mDecs');
    
    %Save 
    assignin('base','no_mObjs',mObjsRet);
    assignin('base','no_mDecs',mDecRet);
    assignin('base','no_vargs',VarsRet);
    ['Figure data saved. To recreate, enter the command:']
    ['nearoptplotmo2(no_mObjs, no_mDecs, no_vargs{:})']
end

function LoadData(hind,evnt,hCurr,dfAns)
    % Allows the user to provide new data in the form of additional
    % objective function, decision variable, and grouping values for
    % additional alternatives to add to the plot.
    
    % INPUTS
    %   hCurr = handle to the current figure window
    %   dfAns = optional cell array of strings with default answers
    %   (default is empty)
    
    [VarsRet, hCont, mObjsOrig, mDecsOrig] = ReadControls('LoadData',hCurr);
    [vGroupOrig,mGroupDataOrig,nO] = aGetField(VarsRet, {'vGroup' 'mGroupData' 'nO'});  
    nD = size(mDecsOrig,2);
       
    %Step 1. Solicit the user for additional alternatives to add
    sDirections = ['Enter data for one or more new alternatives (nNew). Rows in the entries below represent individual alternatives. All entries must be valid Matlab commands and can reference data, functions, or existing variables/matrices in the Base Workspace (note: ', ...
        num2str(nD),' = number of decision variables)', char(10), char(10)];
    sObjs = ['New objective function values (nNew x ',num2str(nO),'):'];
    sDecs = ['New decision variable values (nNew x ',num2str(nD),'):'];
    sGroups = ['Group name(s) [(nNew x 1 cell array of string) | single string (all same group)]:'];
    
    if nO==0
        lNumFields = 2;
        cPrompts = {[sDirections, sDecs],sGroups};
    else
        lNumFields = 3;
        cPrompts = {[sDirections, sObjs],sDecs,sGroups};
    end
        
    options.interpretor = 'tex';
       
    %Fields are:
    %   1 - matrix of new Objective function values (nNew x nO) ***Optional
    %   2 - matrix of new Decision variable values (nNew x nD)
    %   3 - Groupings [cell array of string (nNew x 1) | single string | empty]
       
    if (nargin <= 3) || (length(dfAns) ~= lNumFields)
        dfAns = repmat({''},lNumFields,1);
    end
    drawnow; pause(0.05); % Show a clean dialog
    vAnswers = inputdlg(cPrompts,'Load Data',1,dfAns,options);
    
    %Step 2. Check validity of input
     if isempty(vAnswers) 
        %Canceled, exit immediately
        return
     elseif all(strcmpi(vAnswers,''))
        %No input provided, warn
        h = errordlg('No imput provided on new alternatives. Canceling','Input Error');
        uiwait(h);
        return
    end
    
    %Go though each input and check validity
    cAnswers = cell(lNumFields,1);
    mSizes = zeros(lNumFields,2);
    for i=1:lNumFields
       try 
          cAnswers{i} = evalin('base',vAnswers{i}); 
       catch
          h = errordlg(['Error with input: ',vAnswers{i},'. Try again.'],'Input Error');
          uiwait(h);
          vAnswers{i} = '';
          LoadData(0,0,hCurr,vAnswers);
          return
       end
       mSizes(i,:) = size(cAnswers{i});
    end
    
    %Check the type of the last grouping field
    if mSizes(lNumFields,1) == 1
        %Single entry, expand to all alternatives
        [sNewGroupName, mGroupDataNew] = GenerateNewGroup(hCurr,cAnswers{lNumFields}); 
        cAnswers{lNumFields} = repmat({sNewGroupName},mSizes(1,1),1);
        mSizes(lNumFields,:) = size(cAnswers{lNumFields});
    else
        %List of entires, extract uniques from added groups and insert
        %after first checked group
        
        vUniques = unique(cAnswers{lNumFields});
        lNumNewGroups = length(vUniques);
        cbGroupChecks = aGetField(hCont,{'GroupChecks'});
        [sFirst,lFirst] = GetFirst(cbGroupChecks,mGroupDataOrig(:,1));
        lNumOrigGroups = size(mGroupDataOrig,1);      
        
        mGroupDataNew = [mGroupDataOrig num2cell(ones(lNumOrigGroups,1).*[1:lNumOrigGroups]'); vUniques num2cell([ones(lNumNewGroups,2) lFirst+0.01*[1:lNumNewGroups]'])];
        mGroupDataNew = sortrows(mGroupDataNew,4);
    end
        
    %Check for consistency with existing model inputs
    if nO==0
        lTrueSize = [nD mSizes(1,1)];
        lFields = [1 2];
        iSzCol = [2 1];
        mObjsAdd = [];
    else
        lTrueSize = [nO nD mSizes(1,1) mSizes(1,1)];
        lFields = [1 2 2 3];
        iSzCol = [2 2 1 1];       
        mObjsAdd = cAnswers{1};
    end
    
    for i=1:length(lFields)
        if (lTrueSize(i) ~= mSizes(lFields(i),iSzCol(i)))
            sDir = 'rows';
            if iSzCol(i) == 2
                sDir = 'columns';
            end
            
            h = errordlg(['Error with input: ',vAnswers{lFields(i)},'. Number of ', sDir, ' must be ',num2str(lTrueSize(i)),'. Try again.'],'Input Error');
            uiwait(h);
        	vAnswers{lFields(i)} = '';
            LoadData(0,0,hCurr,vAnswers);
            return
        end
    end      
        
    %Step 3. Update the model formulation
    mObjsNew = [mObjsOrig; mObjsAdd];
    mDecsNew = [mDecsOrig; cAnswers{lNumFields-1}];
    vGroupNew = [vGroupOrig; cAnswers{lNumFields}];

    VarsNew = aSetField(VarsRet,'vGroup',vGroupNew,'mGroupData',mGroupDataNew(:,1:3));       
    
    %Step 4. Replot
    nearoptplotmo2(mObjsNew,mDecsNew,VarsNew{:});
end

function LoadPareto(hind,event,mMatrix,nObjs,txtTolerance,cbChecks,sSliders,mConvert,ckAdd,txtAllowableDeviation,varargin)
    %Loads and replots pareto optimal solutions
    %If ckAdd is checked, adds pareto solutions on top of the existing
    %solutions in mMatrix
    %If ckAdd is unchecked, only plots pareto solutions
    %graph
    
    [mIn nIn] = size(mMatrix);
        
    %Collect the values from the controls
    ToleranceValue = str2num(get(txtTolerance,'String'));
    AllowableDeviationNew = GetCheckTxtValue('LoadPareto',txtAllowableDeviation);
    
    vFixed = ReadCheckboxValues(cbChecks);
    blAdd = get(ckAdd,'Value');
    [m n] = size(vFixed);
    vFixedActions = vFixed((nObjs+1):n);
    
    [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    [vXLabelsShort iOpt vGroup] = aGetField(varargin,{'vXLabelsShort' 'iOpt' 'vGroup'});
    %determine which vGroup is the orginal/base group
    [vUniques, iC, uIndex, iCBeg, iVBeg] = GroupOrderAppear(vGroup);
    lNewGroup = length(vUniques)+1;
    
    %Load in the Pareto optimal results spit out by GAMS to Excel
    mPareto = xlsread('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\UMYMVOut.xls','ToMatlab','AA29:AT46');
    [vNum vLabelsPareto vRaw] = xlsread('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\UMYMVOut.xls','ToMatlab','AC28:AT28');
  
    [mP nP] = size(mPareto);
    
    if nP~=nIn
        fprintf('Current and Pareto data matrixes are different sizes: %d %d\n',nIn,nP);
        return
    end
    
    %convert from JD to $US
    mPareto(:,1) = 1.41*mPareto(:,1);
    mPareto(:,2) = 1.41^2*mPareto(:,2);
    
    nXL = max(size(vXLabelsShort));
    %re-order columns so they are the same as previously plotted data
    [vXSort vISort] = sort(vXLabelsShort);
    [vXSortPareto vISortPareto] = sort(vLabelsPareto);   
    %now move the columns to the appropriate places so they are
    %compatible with mData
    mNewPareto = mPareto;
    for i=1:nXL
        mNewPareto(:,vISort(i)+nObjs) = mPareto(:,vISortPareto(i)+nObjs);
    end

    [mObjsOrig mDecsOrig] = SplitMatrix(mMatrix,nObjs);
    
    mNewParetoNE = mNewPareto;
    
    mRowsToImport = 0;
    %determine number of pareto rows to import, i.e. up to when 
    %   Pareto objective function value <= (Near Optimal Tolerance)*(Optimal
    %   Objective function value)
    
    if iOpt<1
        mRowsToImport = mP;
    else
        for i=1:mP
            if mNewPareto(i,1)<=ToleranceValue*mObjsOrig(iOpt,1)
                mRowsToImport=mRowsToImport+1;
                mNewParetoNE(mRowsToImport,:) = mNewPareto(i,:);
            end
        end        
    end

    [mObjsPar mDecsPar] = SplitMatrix(mNewParetoNE,nObjs);
    
    if blAdd == 1
        %add to original matrix
        mObjs = [mObjsOrig; mObjsPar(1:mRowsToImport,:)];
        mNewMatrix = [mDecsOrig; mDecsPar(1:mRowsToImport,:)];
        %vGroup = ones(mIn+mRowsToImport,1);
        
        if isnumeric(vGroup)
            vGroupNew = [vGroup; lNewGroup*ones(mRowsToImport,1)]; % vGroup((mIn+1):(mIn+mRowsToImport),1) = 2;  
        else
            vGroupNew = vGroup;
            for i=1:mRowsToImport
                vGroupNew{mIn+i} = sprintf('%d. Pareto optimal',lNewGroup);
            end
        end
    else
        %Only show pareto optimals
        mObjs = mObjsPar(1:mRowsToImport,:);
        mNewMatrix = mDecsPar(1:mRowsToImport,:);
        if isnumeric(vGroup)
            vGroupNew = ones(mRowsToImport,1);  
        else
            for i=1:mRowsToImport
                vGroupNew{i} = sprintf('%d. Pareto optimal',1);
            end
        end
         
        iOpt=0;
    end
    
    fprintf('%.0f Pareto solutions added\n',mRowsToImport);
    vGroupNew;
    mObjsOrig;
    mObjsPar;
    
    mObjs;
    mNewMatrix;
    vFixedVals;
    varargout = aSetField(varargin,'tolerance',ToleranceValue,'vFixed',vFixed,'vFixedVals',vFixedVals,'vStep',vStep, 'GroupToHighlight',iOpt,'vGroup',vGroupNew, 'SubSapceError', AllowableDeviationNew);
    nearoptplotmo2(hWind, mObjs, mNewMatrix, varargout{:});
end

function TestAllRegion(hind,event,hWind,mMatrix,nObjs,iRunGAMS)
    %called when user presses the Test All region button
    %Reads in Tolerance and check box values and writes to a text file GAMS
    %can read. Then determines all permutations for free (unchecked) axes
    %and replots into the figure
    
    % INPUTS
    %  hWind = handle to figure object to pull out stored variables
    %  mData = the data in data units
    %  nObjs = number of objectives (columns) in the data set
    %  hControls = cell array of handles to controls on the figure
    %  mConvert = matrix of axis conversion factors from plot to data units
                            
    % Get the handles for the figure controls we will need
    [varargin_fresh,hControls] = ReadControls('TestAllRegion',hWind);
           
    %Read in the variable arguments
    [ToleranceValue,NumSamples,vFixed,vFixedVals,vXLabelsShort,vGroup] = aGetField(varargin_fresh,{'tolerance','NumSamples','vFixed','vFixedVals', 'vXLabelsShort','vGroup'});
    [m n] = size(vFixed);
    
    %Calc number of decision variables
    nD = n-nObjs;
    
    %Collect the values from the controls
    vFixedActions = vFixed((nObjs+1):n);
    
    [mObjs, mData] = SplitMatrix(mMatrix,nObjs);

    vMaxResult = max(mData);
    vMinResult = min(mData);
    
    if iRunGAMS==1
        %Calculate the number of solutions the selected fixed and non-fixed
        %variables and ranges can generate
        dNumSols = ((1-vFixedActions).*((vMaxResult - vMinResult)./vStep+1));
         
        %find last non-zero element
        for i=n-nObjs:-1:1
            if dNumSols(i) > 0
                lLastElement = i;
                break;
            end
        end
        
        %convert zeros to ones
        for i=1:n-nObjs
            if dNumSols(i)==0
                dNumSols(i)=1;
            end
        end
        
        %sort the number of sols into increasing order
        [dSols iSols] = sort(dNumSols);
        
        %calculate the # of +/- branches
        dNumRangeSols = 0;
        for i=1:lLastElement-1
            if dNumSols(i) > 1
                dNumRangeSols = dNumRangeSols + 2*prod(dNumSols(1,1:i));

            
%            if dSols(i) > 1
%                dNumRangeSols = dNumRangeSols + 2*prod(dSols(1,1:i));
            end
        end
        
        %sort the # of 
        dNumSols = prod(dNumSols);      

        %write to the file
        fInt = fopen('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\matlab.inc','w');
        fprintf(fInt,'*  Values set by user in MatLab\n');
        fprintf(fInt,'   SET sol Near Optimal Solutions /sol1*sol%.0f/;\n', dNumSols);
        fprintf(fInt,'   SET rngSol Range solutions /rsol1*rsol%.0f/;\n', dNumRangeSols);
        fprintf(fInt,'   GAMMA = %s;\n',ToleranceValue);
        fprintf(fInt,'   FOPT = %.2f;\n', mData(1,iOpt)/1.41);
        %fprintf(fInt,'PARAMETER LFIXED(i)  Long-term action is fixed (1 = fixed and 0=free) /\n');
        for i=1:n-nObjs
            fprintf(fInt,'   LFIXED(''%s'')= %i;\n',char(vXLabelsShort(i)),vFixedActions(i));
        end
        for i=1:n-nObjs
            fprintf(fInt,'   LFixedValVol(''%s'')= %.3f;\n',char(vXLabelsShort(i)),vFixedVals(i+nObjs));
        end
        for i=1:n-nObjs
            fprintf(fInt,'   XMINs(''%s'')= %.3f;\n',char(vXLabelsShort(i)),vMinResult(i));
            fprintf(fInt,'   XMAXs(''%s'')= %.3f;\n',char(vXLabelsShort(i)),vMaxResult(i));
        end
        for i=1:n-nObjs
            fprintf(fInt,'   cMAP(''%s'')= %.0f;\n',char(vXLabelsShort(i)),iSols(i));
        end

        fclose(fInt);
   
        %code to do iterative permulations of solution values

        vXStart = vFixedVals((nObjs+1):n);

        [m2 n2] = size(vXStart);
        if m2>n2
            vXStart = vXStart';
        end

        mNewResults=ExploreSubregion(1,vXStart,vMinResult,vMaxResult,vFixedActions,vStep,1);

        [nM nN] = size(mNewResults);

        %write to file to read in by gams
        fInt = fopen('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\sol_values.inc','w');
        fprintf(fInt,'*  Solution combinations determined by MatLab based on fixed variable combinations\n');
        fprintf(fInt,'   PARAMETER XSOLs(sol,i) Pre-specified solution combinations (MCM per year);\n');
        for i=1:nM %solutions
            for j=1:nN %decision variable values
                fprintf(fInt,'   XSOLs(''sol%.0f'',''%s'')= %.2f;\n',i,char(vXLabelsShort(j)),mNewResults(i,j));
            end
        end

        %now run GAMS
        [status,result] = dos('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\rumyneoptr.bat &','-echo');
        h = msgbox(sprintf('Running %d solutions\nPress OK once GAMS finishes running',nM),'Running GAMS');
        uiwait(h);
    else
        %read solutions from prevous run       
        mActVolumesRead = dlmread('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\mat_nevol.csv',',',1,0);
        [nM nN] = size(mActVolumesRead);
        mActVolumesRead=mActVolumesRead(:,1:(nN-1));
        nN = nN-1;
        
        %read in column headers for solutions from previous run
        [vLabelsRead vFixedRead] = textread('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\mat_neact.csv','%s%d','headerlines',1','delimiter',',');
        [vXSort, vISort] = sort(vXLabelsShort);
        [vXSortRead vISortRead] = sort(vLabelsRead);
        
        %now move the columns to the appropriate places so they are
        %compatible with mData
        mNewResults = mActVolumesRead;
        vFixed = vFixedRead;
        
        for i=1:nN
            mNewResults(:,vISort(i)) = mActVolumesRead(:,vISortRead(i));
            vFixed(vISort(i)) = vFixedRead(vISortRead(i));
        end        
        
        %transpose vFixed and add zeros for objective functions
        vFixed = [zeros(1,nObjs) vFixed'];
    end
    
    %Load GAMS output back in
    [Year MVI Sol Total Variance ModStat SolStat] = textread('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\mat_neobj.csv','%d%s%s%f%f%d%d','headerlines',1','delimiter',',');
    
    %convert from JD to $US
    Total = 1.41*Total;
    Variance = 1.41^2*Variance;
    
    blNeedFirst = 1;
    %post process results - we only want solutions where SolveStatus = 1 and Model Status = 1, 2, 8
    mFinalResult = mData;
    mFinalObjs = mObjs;
    
    %size(mFinalResult);
    %size(mFinalObjs);
    %size(Total);
    %size(Variance);
    
    for i=1:nM
       if SolStat(i)==1 && ((ModStat(i)==1) || (ModStat(i)==2) || (ModStat(i)==8))
           %keep this solution
           mFinalResult = [mFinalResult; mNewResults(i,:)];
           cObjRow = [Total(i) Variance(i)];
           mFinalObjs = [mFinalObjs; cObjRow(1:nObjs)];
       end
    end

    [mF nF] = size(mFinalResult);
    [mI nI] = size(mData);
    
    %add a new group for the newly loaded values. This group will be added
    %after the first checked group
    if iscell(vGroup)
        [sNewGroup, mGroupDataNewSort] = GenerateNewGroup(hWind,'NearOpt');
        vGroupNew = [vGroup; repmat({sNewGroup},mF-mI,1)];
        %vGroupNew = vGroup;
        %for j=mI+1:mI+mF
        %    vGroupNew{i} = sNewGroup;
        %end
    else
       [sNewGroup, mGroupDataNewSort] = GenerateNewGroup(hWind,'NearOpt'); 
       vGroupNew = [vGroup; sNewGroup*ones(mF-mI,1)];
    end   
    
    [min(mNewResults);max(mNewResults)]
    fprintf('%d feasible near-optimal solutions found and added\n',mF-mI);
        
    %add objective function row to vFixedActions and vFixedValues
    %vFixedValues = [vFixedVals(1); vFixedValues];
    %vActLong = ['f(X)' ; vActLong];
     
    varargout = aSetField(varargin_fresh,'vGroup',vGroupNew,'mGroupData',mGroupDataNewSort(:,1:3));
    nearoptplotmo2(mFinalObjs, mFinalResult, varargout{:});
end

function [blReturn] = Enabled2Boolean(sEnabled,lToggle,cStr)
    %Takes the string input sEnabled (on/off) and returns the corresponding
    %integer value (1/0)
    % set optional lToggle to 1 to return the opposite of the incoming
    % value (Default value 0 [no toggle])
    % cStr is an option cell array of strings with the categories to test. Default value is {'off' 'on'}. 
    
    if nargin < 2
        lToggle = 0; %no toggle
    end
    if nargin < 3
        cStr = {'off' 'on'}';
    end
    
    %set the toggle mapping
    if lToggle ==1
        a = -1;
        b = 1;
    else
        a = 1;
        b = 0;
    end
    
    blReturn = a*(find(strcmpi(cStr,sEnabled)==1,1,'first')-1)+b;
end

function [sReturn] = Boolean2Enabled(blValue,lToggle,cStr)
    %Takes the boolean or integer input blValue and returns the corresponding
    %integer string value
    %set optional lToggle to 1 to return the opposite of the incoming
    % value (Default value 0 [no toggle])
    % cStr is an option cell array of strings with the categories to return. Default value is {'off' 'on'}. 
    if nargin < 2
        lToggle = 0; %no toggle
    end
    if nargin < 3
        cStr = {'off' 'on'}';
    end
    
    %set the toggle mapping
    if lToggle == 1
        a = -1;
        b = 1;
    else
        a = 1;
        b = 0;
    end
    
    sReturn = cStr{a*(blValue+1)+b};
end

    