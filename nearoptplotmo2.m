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
%       A row in each matrix mObjs and mDecisions represents a solution to the underlying optimization
%       problem (i.e., a "point" in the multivariate space).
%
%   2) Tools on the 'Display' tab (far right) to brush, pivot, and interactively manipulate the objective
%       function and decision variable data matricies as well as rows within them (i.e., both by columns/axes and rows/groups/color).
%
%   3) Tools on the 'Interact' tab to generate new solutions or
%       alterantives, e.g., for near-optimal analysis. There are three data sources available to use to generate
%       new solutions: i) the data matrix itself (querying), 2)
%       MATLAB (for linear programs using data specified in the objective
%       function vector(s) and constraint matrix), and 3) executing a GAMS
%       file (note, you will need to install GAMS, see www.gams.com). The
%       methods can generate: a) one solution, b) the maximum extents
%       representing the minimum and maximum values given current selected
%       variables and their values, c) randomly generated solutions, and d)
%       all solutions via enumeration (for MIP problems)
%
%   4) Tools on the 'File' tab to save the current figure settings to the Matlab base workspace and use
%       to recreate the figure
% 
%   5) A variety of parameters/settings to control how the parallel
%       coordinate plot is displayed and labeled (varargin)
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
%      and is followed by the parameter value.
%      If a parameter is omitted or in incorrect format, the default value
%      is used. e.g., varargin = {'hWind',1,'SubSpaceError',0.2} would set the
%      parameters hWind and SubSpaceError, to the values, respectively,
%      of 1 and 0.2. Below is a listing of the parameters followed
%      by their meanings and default values (listed in parenthesis).
%
%     hWind = handle for figure window in which to plot. (Default: 0 = plots in a
%         new figure)
%
%     tolerance = near optimal tolernace. Fraction of the optimal objective function value and used to generate
%         new alternatives for optimization problems. Populates a text box
%         on the Interact tab for generating additional data sets. (Default
%         value is 1)
%
%     SubSpaceError = the error allowed in selecting data points within a
%         (small) range of the set value to highlight on the current plot
%         during mouseover or include on a sub-space plot. This value is specified in
%         units of fraction of a tick and can
%         also be dynamically changed by entering a value in the textbox on
%         the Interact Tab. (Default is 0 [no error]).
%         E.g., if the ticks for the axis are in units of 0, 5000, 10000,
%         ..., entering a SubSpaceError of 0.1 will select all solutions
%         within +/- 500 of the set value. This value is also used when
%         re-sampling new points within the polytope defined by AMat and
%         Bhrs.
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

%     fontsize = fontsize to display axis labels and all associated
%            deravative labels, controls, etc. (Default value: 16)
% 
%     mActCat = an nD by nF matrix whose rows represent the decision variables and nF columns are grouping fields to categorically classify
%           the decision variables; used to generate an optional table of decision variable labels to position below the decision variable axis labels and to
%           color-code groups of those labels/axes (default value - [])

%     vGroup  = m-length column vector whose values specify the group to
%           which each row in mObjs and mDecisions belongs and should be
%           colored on the plot. Can either be a numeric or text string. The unique values of vGroup are used to generate
%           the group labels and checkboxes on the Display tab where the user can dynamically specify which groups to show on the plot
%           (Default value: all ones [1 group])
%
%     vGroupOrder = an nU length column vector whose entires (values in vGroup) specify the
%           order in which to plot groups. The first group listed is plotted as the lower-most
%           layer and will appear first in the group list on the Display
%           tab (Default value: alpha-numeric sort)
%
%     mGroupData = an nU by 2 cell matrix whose entries describe attributes of each unique group in vGroup. 
%         Column 1 - vGroupName - The name of the group. The first group listed is plotted as the lower-most
%                   layer and will appear first in the group list on the
%                   Display Tab.
%         Column 2 - ShowGroup - binary variable (1=yes) to show the group
%                   on the plot and have the checkbox next to the group
%                   name checked.
%       (Default values -- alphanumeric sort and all groups checked)
%    
%     GroupToHighlight = value in vGroup representing a group to highlight on
%           the plot with thick lines of the color mHighlighColor (see below). E.g., use this feature to 
%           highlight the optimal solution (Default: 0 [no group])

%     CurrentRecord = row in mObjs and mDecisions to highlight on the plot.
%           (Default: 1 or the record number of a record in
%           GroupToHighlight if a GroupToHighlight is specified)%     
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
%     mHighlightColor = a 1 by 3 matrix indicating the colors to use for
%           highlighted lines on the plot (Group values =
%           GroupToHightlight). (Default value is Thick Black for 
%           all [0 0 0])
%
%     SubSpaceGroup = value in vGroup representing the current group of solutions in the highlighted subspace (Default: 0 [none])
    
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
%     AMat = lCon by nD matrix of coeffiecients for the lCon constraints
%           in the underlying linear optimization problem associated with the
%           data provided. This parameter, along with Brhs and possibly cFunc, is used when clicking the 'Random Sample' button on the
%           sliders tab to sample new random solutions (default value [])
%
%     Brhs = column vector of lCon elements whose values are the right hand side
%           of the constraints that define the system of equations from which to randomly sample additional points.
%           i.e., AMat * X <= Brhs. (Default value [])
%
%     cFunc = nO by nD matrix of coeffiecents where each row respresents
%           the coefficients for the multiple objective functions of the
%           optimization problem. Used with AMat and Brhs to sample
%           subspaces of the solution set. (Default value [])

%     ShowControls = boolean takes the value of 1 to show the plot
%           controls as menu items on the right and checkboxes/sliders
%           (Default Value: 1 [show controls])
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
%           (Default of [0.100 0.615 0.47 0.6182 - yBottom + 0.3068])
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
% D) nearoptplotmo2(150000*rand(15,1),[5*rand(15,1) 20*rand(15,1) 5*rand(15,2)],'vObjLabels',{'Cost ($)'},'vXLabels',{'Area A' 'Area B' 'Machine 1' 'Machine B'},'AxisScales','auto');
%       Labels the objective function and decision axes. Default AxisScales
%       setting gives similar scaling as Case C.
%
% E) nearoptplotmo2(150000*rand(30,2),[5*rand(30,1) 20*rand(30,1) 5*rand(30,6)],'vGroup',[ones(20,1);2*ones(10,1)],'AxisScales','auto');
%       Displays 30 traces across 9 axes. Assigns the first 20 rows to Group 1 plotted in light green and
%       last 10 rows to Group 2 plotted in blue.
%
% F) MyGroups = {'Group 2' 'Group 2' 'Group 2' 'Group 1' 'Group 1'}';
%    nearoptplotmo2(150000*rand(5,2),[5*rand(5,1) 20*rand(5,1) 5*rand(5,6)],'vGroup',MyGroups,'AxisScales','auto');
%       Like E but text labels used for groups. Note groups are ordered
%       alphanumerically so Group 1 still listed and plotted first in light green.
%   
%
%   Programmed by David E. Rosenberg
%   July 2012
%   Updated May 2013 to include passing optional input arguements as a variable list and
%      using the sliders and axes checkboxes to view subsets of solutions
%   Updated June 2014 to include solution generation and interaction functions, reorganize GUI, allow data save, 
%      and improve documentation.
%   
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   This code and all the dependent functions are distributed as freeware AS-IS with
%   no expressed or implied warrenty regarding the claimed functionility. To report bugs,
%   email david.rosenberg@usu.edu (note, there is no promise of when, if
%   ever, the bug will be corrected).

%% Read in stuff about the main input parameters
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
    Tolerance = 1;
    SubSpaceError = 0;
    NumSamples=0;
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
         vXLabels{i} = ['Decision ', num2str(i)];
    end
    
    vXLabelsShort = vXLabels;
    yAxisLabels = {'Objective Functions' 'Decision Variables'};

    fontsize = 16;
    mActCat = ones(nD,1);
    AxisScales='none';
    BaseAxis=[min(nO,1) 1];  %put a zero in the first element if there are no objective axes
    NumTicks = [6 6];
    TickMag = 5;
    mLims = zeros(2,n);
    
    iOpt = '0';
    lCurrRecord = 1;
    SubSpaceGroup = 0;
    
    vGroup = ones(m,1);
    vGroupOrder = [];
    mGroupData = {};
    mColors = [];
    mHighlightColor=[0 0 0];
    
    vMaxValues = max(mDecisions);
    AMat = []; AMatTemp=[];
    Brhs = []; BrhsTemp=[];   
    cFunc = []; cFuncTemp=[];
    
    HideSliders = 1;
    GenType = 1;
    GenMethod = 1;
    sGamsFile = '';
    
    
    yBottom = 0.47;
    PlotPosition = [0.100 yBottom 0.615 0.6182 - yBottom + 0.3068];
    PanelWidth = 285;
    
    ShowControls = 1;
    
    ButtonText = {'File' 'Interact' 'Display'};
    StartTab = 3;
    lStartTab = StartTab; %Integer value
    cStartTab = ButtonText(lStartTab); %Cell value
    
    %Fields to ignore
    IgnoreFields = {'mConvert' 'SubSpaceErrorPU'};
    
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
            
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'tolerance'))            
                Tolerance = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'SubSpaceError'))            
                SubSpaceError = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'NumSamples'))            
                NumSamples = varargin{count+1};
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
                varargin{count};
                varargin{count+1};
                
                if (size(varargin{count+1},1) == nD)
                    mActCat = varargin{count+1};
                else
                    warning(['nearoptplotmo2: mActCat has a different number of rows than the number of columns in mDecisions', num2str(nD), '. Continuing with default mActCat.'])
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
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'fontsize'))
                fontsize = abs(varargin{count+1});
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
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mColors'))
                mColors = varargin{count+1};
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'mHighlightColor'))
                if (size(varargin{count+1},2) == 3)
                   mHighlightColor = varargin{count+1};
                else
                    warning(['nearoptplotmo2: mHighlightColor must only be a 1 by 3 matrix. Continuing with default highlight color setting.'])
                end
                
                count=count+1;
                
            elseif (ischar(varargin{count}) && strcmpi(varargin{count},'SubSpaceGroup'))
                SubSpaceGroup = varargin{count+1};
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
                    warnmsg = ['nearoptplotmo2: cFunc is improperly sized. Default to empty [].'];
                    cFuncTemp = [];
                    warning(warnmsg)
                end
                count=count+1;
                
         
             elseif  (ischar(varargin{count}) && strcmpi(varargin{count},'PlotPosition'))
                PlotPosition = varargin{count+1};
                
                if (length(varargin{count+1})==4)
                    PlotPosition = varargin{count+1};
                else
                    warnmsg = ['nearoptplotmo2: PlotPosition has', num2str(length(PlotPosition)), ' but needs 4. Continuing with default position.'];
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
                                   
                    if (exist([varargin{count+1},'.gms'],'file'))
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
        if (length(BaseAxis) ~= 2) || (BaseAxis(1)>nO) || (BaseAxis(2)>nD)
            warning(['nearoptplotmo2: AxisScale BaseAxis does not have 2 elements or values are larger than the number of axes. Continuing with default values of 1st objective and decision axes.'])
            BaseAxis = ones(1,2);
        end
        if (size(mLims,2) ~= nO+nD) || (size(mLims,1)<2)
            warning(['nearoptplotmo2: AxisScale mLims is improperly sized. Continuing with default AxisScales setting of auto.'])
            AxisScales='auto';
        end
        count=count+3;
        
        %Further error checking on AMat and Brhs          
            {class(AMatTemp) class(BrhsTemp)};
            AMatTemp;
            BrhsTemp;
           
            
            if ~isempty(AMatTemp) || ~isempty(BrhsTemp)
                if ~isempty(AMatTemp)
                   if (size(AMatTemp,2) == nD)
                        if (size(AMatTemp,1) == size(BrhsTemp,1))
                             AMat = AMatTemp;
                             Brhs = BrhsTemp;
                             cFunc = cFuncTemp;
                        else
                             warning(['nearoptplotmo2: # rows of AMat and Brhs are different.', ...
                        ' Continuing with default of empty inputs.'])
                        end
                    else
                        warning('nearoptplotmo2: # columns of AMat must equal number of columns of mDecision. Reverting to default empty AMat, Brhs, cFunc inputs')
                    end
                 else
                    warning('nearoptplotmo2: need to define both an AMat and Brhs. No AMat. Reverting to default empty inputs for both and cFunc.')
                end
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
          lHighRecord = lCurrRecord; %Also save as the highlighted record
       else
          if lCurrRecord < 1
              lCurrRecord = 1;
          elseif lCurrRecord > m
              lCurrRecord = m;
          end             
       end
    
    
        
    %Determine the number of groups.
    [vUnique vIC vID] = unique(vGroup);
    nU = max(size(vUnique));
    vShowGroup = ones(nU,1);
    
    if ~isempty(mGroupData)
        
       %vUnique;
       %mGroupData(:,1);
        
       %Compare vUnique to vGroupOrder 
       IsMemberInd = ismember(mGroupData(:,1),vUnique);
       if (sum(IsMemberInd)==length(vUnique)) && (length(vUnique)==size(mGroupData,1))
           vUnique = mGroupData(:,1);
           vShowGroup = (cell2mat(mGroupData(:,2))>0);
       else
           warning('mGroupData has different members than vGroup. Reverting to default alpha-numeric sorting of groups and showing all groups')
       end
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

    mColors;
    
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
    
    %mColors;
    
    % Error checking on the precision -- round to nearest integer
    Precision = round(Precision);
    
    % Round the mResults to the specified precision
    mResults = round(mResults*10^Precision)/10^Precision;
    
    % update vStep if no value passed
    if isnan(vStep)
        vStep = (max(mResults)-min(mResults))/20;
    end
    
    %% Start building the plot
     
    %Get the Figure Handle
    if (hWind>0) && ishandle(hWind)
        Figure2 = hWind;
        clf(Figure2);
    else
        Figure2=figure; %('units','normalized','outerposition',[0 0 1 1]);
        %set(Figure2,'units','pixels');
    end
    
    hWindReturn = Figure2;
       
    lLineWidth = 1.5;
    
    %Determine the size of the figure components (Parallel Plot and Control
    %Panel)
    %Shift to pixel units to get current figure size
    cUnits = get(hWindReturn,'Units');
    set(hWindReturn,'Units','pixel');
    vSizePx = get(hWindReturn,'Position');
    set(hWindReturn,'Units',cUnits);
    
    lMarginPanel = 0.095; % margin between parallel coordinate plot and control panel in normalized units
    rMarginPanel = 0.01; % margin between control panel and figure right edge in normalized units
    xPanelWidthNorm = PanelWidth/vSizePx(3); %width of control panel in normalized units
    
    %Covert the plot position to local variables
    xLeft = PlotPosition(1);
    xWidth = 1 - (xLeft + lMarginPanel + xPanelWidthNorm + rMarginPanel); %PlotPosition(3);
    yBottom = PlotPosition(2);
    yHeight = PlotPosition(4);
    
    
    
    %xLeft = 0.100;
    %xWidth = 0.615;
    %yBottom = 0.47; %0.3068;
    %yHeight = 0.6182 - yBottom + 0.3068; %0.6182;
    
    %plot2 = subplot(1,1,1,'FontSize',fontsize-4,'Position',[0.13 0.2429 0.6593 0.6821]);
    
    if nO==0
        vColorLeftScale = mColors(1,1,2,:);
    else
        vColorLeftScale = mColors(1,2,2,:);
    end
    
    plot2 = axes('Parent',Figure2,'FontSize',fontsize-4,'Units','normalized','Position',[xLeft yBottom xWidth yHeight],'YColor',vColorLeftScale,'box','on');
    
    hold on
       
    %fprintf('Initial Maxes / Mins\n');
    vMaxResult = max(mResults);
    vMinResult = min(mResults);
    
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
    
    %set the limits and get the ticks   
    if 1
        set(plot2,'xLim',[1 n+0.5]);
        %dTemp = plot([ymin ymaxdec]);
        %yLimsStart = get(plot2,'yLim');
        %dticks = [yLimsStart(1):(yLimsStart(2)-yLimsStart(1))/(NumTicks-1):yLimsStart(2)];
        %delete(dTemp);
    else %old approach
        dticks = [ymin:(ymaxdec-ymin)/(NumTicks(1)-1):ymaxdec];
        set(plot2,'yLim',[ymin ymaxdec],'xLim',[1 n+0.5],'ytick',dticks,'yticklabel',dticks);
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
                %calc for just the objective functions
                %[mUseResObj, vMultObj, BaseObjCol] = ScaleMultipleAxes(mPlotTemp(:,1:nO));
                
                %MaxObjVal = max(mUseResObj(:,BaseObjCol));
                %BaseAxis(1) = BaseObjCol;
                %BaseObjCol;
                %vMultObj;
                %Now build the entire set of transformed objectives and decision axes
                [mUseResAll, vMultAll, BaseAll] = ScaleMultipleAxes(mPlotTemp);
                
                mPlot = mUseResAll;
                BaseAll;
                
                MaxAllVal = max(mUseResAll(:,BaseAll));
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
        if strcmpi(AxisScales,'auto')
            %currticks;
            minDC;
            %YMaxPU = ConvertPlot2DataUnits(max(rticksnum),mTransformToOrig(:,nO+minDC),1);
            %YMaxPU = MaxToUse;
            %YMinPU = MinToUse;
        else
            %YMaxPU = ymaxdec;
            %YMinPU = ymin;
        end
    end
    mTransformToOrig;    
    
    %Revisit the scaling now in plot units   
    if 1
        %Temportary plot in plot units to get the auto Y limits
        if 0
            [min(mResults)' max(mResults)']
            [min(mPlot)' max(mPlot)']
            min(min(mResults))
            max(max(mResults))
        end 
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
    else %old approach
        set(plot2,'xLim',[1 n+0.5],'yLim',[YMinPU YMaxPU]);

        get(plot2,'ytick');
        get(plot2,'yticklabel');

        %ymaxobj = max(max(mPlot(:,1:nO))); 
        ymax = YMaxPU;
        ymin = YMinPU;  
    end
    
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
                         [s1 s2] = strread(vAxisLabelsAll{i},'%s %s','delimiter','(');
                         s2 = [' (',s2];
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
        vMaxPlot = max(mPlot);
        vMinPlot = min(mPlot);
   % end
    
    mPlot(:,1:nO);
    %fprintf('Final Maxes/Mins\n');
    vMaxResult;
    vMinResult;
    
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
    
    %mColorChoices = [0.142 0 0.85; 0.65 0 0.13;0 0.737 0]; % blue-red-green
        
    mColorMatrix = zeros(5,5,3);
    mColorMatrix(1,:,:) = [0.6	0.06	0.06; 0.7	0.175	0.175; 0.8	0.32	0.32; 0.9	0.495	0.495; 1	0.7	0.7]; %red
    mColorMatrix(2,:,:) = [0.6	0.33	0.06; 0.7	0.438	0.175; 0.8	0.56	0.32; 0.9	0.697	0.495; 1	0.85	0.7]; % orange
    mColorMatrix(3,:,:) = [0.42	0.6	0.06; 0.525	0.7	0.175; 0.64	0.8	0.32; 0.765	0.9	0.495; 0.9	1	0.7]; % green
    mColorMatrix(4,:,:) = [0.06	0.42	0.6; 0.175	0.525	0.7; 0.32	0.64	0.8; 0.495	0.765	0.9; 0.7	0.9	1]; % blue
    mColorMatrix(5,:,:) = [0.15	0.06	0.6; 0.262	0.175	0.7; 0.4	0.32	0.8; 0.562	0.495	0.9;0.75	0.7	1]; %purple

    %instead build orange-purple-red
    mColorChoices = [mColorMatrix(2,3,:);mColorMatrix(5,3,:);mColorMatrix(1,3,:)];
    mColorFull = [mColorMatrix(2,:,:);mColorMatrix(5,:,:);mColorMatrix(1,:,:)];
    
    mActCat;
    
    [vUniques, iC, uIndex, iCBeg, iVBeg] = GroupOrderAppear(mActCat(:,1));
        
    lMainClass = 1;
    mTextPos = zeros(n,4);
    
    %display the axis lines and labels
    for i=1:n    
           ['i: ', num2str(i),'; nO: ', num2str(nO)];
           
            h3 = plot([i i],[ymin ymax]);
            if i>nO %decisions
                lLineColor = [0.8 0.8 0.8];
                lMainClass = iVBeg(i-nO);
                vColor = mColorChoices(iVBeg(i-nO),:);
                %hTexts(i) = text(i,ymin-6*(ymax-ymin)/100,vXLabels(i-nO),'fontsize',fontsize-6,'HorizontalAlignment', 'right','color',vColor,'VerticalAlignment','middle','Rotation',90);
            else
                lLineColor = [0.6 0.6 0.6];
                vColor = [.737 0 0.737];
            end
            set(h3,'linewidth',0.5,'color',lLineColor);
            hTexts(i) = text(i,ymin-6*(ymax-ymin)/100,vAxisLabelsAll{i},'fontsize',fontsize-6,'HorizontalAlignment', 'right','color',vColor,'VerticalAlignment','middle','Rotation',90);
           
            strCurrUnits = get(hTexts(i),'Units');
            set(hTexts(i),'Units','normalized');
            mTextPos(i,:) = get(hTexts(i),'Extent');
            set(hTexts(i),'Units',strCurrUnits);
    end
    
     
     %mTextPos
    %Alternative way of labeling decision variable axes based on mActCat
    %old data units
    %[hTextGroup, rDex, hGLines] = GenerateGroupLabels(mActCat,1,0,[1  ymin-6*(ymax-ymin)/100 3.3*(ymax-ymin)/100 0 4*(ymax-ymin)/100],'data',fontsize-2);
    
    yBottomTable = 0.025;
    mMaxExtent = min(mTextPos(:,2));
    yTableHeight = max([yBottom-abs(yHeight*mMaxExtent)-yBottomTable-0.005; 0.02*(size(mActCat,2)-1); 0.001]);
    get(plot2,'Position');
    %[yBottomTable mMaxExtent yTableHeight]
    %[hTabComponent hTabContainer]=GenerateGroupLabelsTable(hWindReturn,mActCat,3,[xLeft+(nO-1+0.5-.0025)/(n-0.5)*xWidth yBottomTable xWidth*(n-nO+0.25)/(n-0.5) yBottom-yBottomTable-0.035],'normalized',fontsize);
    %position the table below the largest axis text label
    %[xLeft+(nO-1+0.5-.0025)/(n-0.5)*xWidth yBottomTable xWidth*(n-nO+0.25)/(n-0.5) yTableHeight];
    [hTabComponent hTabContainer]=GenerateGroupLabelsTable(hWindReturn,mActCat,3,[xLeft+(nO-1+0.5-.0025)/(n-0.5)*xWidth yBottomTable xWidth*(n-nO+0.25)/(n-0.5) yTableHeight],'normalized',mColorFull,fontsize-5);

    hTextGroup = []; rDex = []; hGLines = [];
    
    
    [nActR nActC] = size(mActCat);
    yFirstHeight = 0.25*yBottom;
  % [hTextGroup, rDex, hGLines] = GenerateGroupLabels(mActCat,1,0,[xLeft+nO/(n+0.5)*xWidth yBottom-.05 (yBottom-yFirstHeight)/(nActC-1) yFirstHeight 0.025 xWidth],'normalized',fontsize-2);
    
    hTGSize = max(size(hTextGroup));
    for i=1:hTGSize
       cPos = get(hTextGroup(i),'Position');
       %sprintf('%d %s %d %.2f\n',i,get(hTextGroup(i),'string'),cPos(1),cPos(2)) 
    end
    
    blHasPareto = 0;
    
    %nOpts = max(size(iOpt));
    
    %plotting

    lOptGroupCurr = 0;
    
    if 1 %nU > 1
        %overplot additional groups in a separate color/formatting
        %first pull out the group
        
        sprintf('Plotting groups\n');

        %[vSort vSortI] = sort(vGroup); 
        %if isnumeric(vGroup)
        %    vGroupCount = histc(vGroup,vUnique); 
        %else
        %    [uniques vGroupCount] = count_unique(vGroup);
        %end
        
        % vGroupCount = histc(vIC,vID); 
        
%        if (nU == vUnique(nU)) && (vGroupCount(2) > 0)
%            % Has Pareto group
%            blHasPareto = 1;
%        end
       

        %vStart = 1;
        %mColorsDecsGroup = [GreenToMagentaRamp(5,:); 0 0 0.8; 0.8 0 0; 0 0 0]; %[0.526 1 0.526; 0 0.737 0;0 0.316 0]
        %mColorsObjsGroup = [GreenToMagentaRamp(11,:);0	0 0.8; 0.8 0 0; 0 0 0]; %[1 0.526 1;0.947	0 0.947; 0.526 0 0.526];
        
        %size(mColors)
        
        mColorsDecsGroup = mColors(:,1,1,:);
        mColorsObjsGroup = mColors(:,1+(nO>0),1,:);
        
        vLineStyle = {'-' '-' '-' '-' '-' '-' '-'};
        %vLineStyle = {'-' '-' '-.'};
        %hs = zeros(1,2*nU);
        
        %fprintf('Plotting groups\nGroup\tStart\tEnd\n');
        
        hPCGroup = cell(nU,2);
        
        
        for i=1:nU
            %fprintf('%.0f\t%.0f\t%0.f\n',i,vStart,vStart+vGroupCount(i)-1);
            
            if 0
            for j=vStart:vStart+vGroupCount(i)-1
                if j==vStart
                    mGroup = mPlot(vSortI(j),:);
                else
                    mGroup = [mGroup; mPlot(vSortI(j),:)];
                end
            end
            end
            
            [i size(vUnique) size(vGroup) size(mPlot)];
            
            if iscell(vUnique(i))
                blOptGroup = strcmpi(vUnique(i),iOpt);
                mGroup = mPlot(strcmpi(vGroup,vUnique(i)),:);
            else
                blOptGroup = vUnique(i) == iOpt;
                mGroup = mPlot(vGroup==vUnique(i),:);
            end
            
            if i==2
                mGroupSave=mGroup;
            end
            
            lIndToUse = 1+mod(i-1,size(mColors,1)); %vUnique(i)
            
            if blOptGroup
                lOptGroupCurr = i;
                lLineWidth=2;
                mColorsDecs = mHighlightColor;
                mColorsObjs = mHighlightColor;
            else
                lLineWidth=1;
                mColorsDecs = mColorsDecsGroup(lIndToUse,:);
                mColorsObjs = mColorsObjsGroup(lIndToUse,:);

            end
            
            if vShowGroup(i) > 0
                strVis = 'on';
            else
                strVis = 'off';
            end
                       
            hPCGroup{i,1}= parallelcoords(mGroup,'color',mColorsDecs,'LineStyle',vLineStyle{lIndToUse},'linewidth',lLineWidth,'Standardize',UseStandarize,'Visible',strVis);
            if nO>0
                hPCGroup{i,2}= parallelcoords(mGroup(:,1:(nO+1)),'color',mColorsObjs,'LineStyle',vLineStyle{lIndToUse},'linewidth',lLineWidth,'Standardize',UseStandarize,'Visible',strVis);
            end
            %vStart = vStart+vGroupCount(i);          
        end
    else
        %no grouping -- plot all rows in the same color
        
        lLineWidth=1;
        %vGroup = ones(1,m);
        h2 = parallelcoords(mPlot,'color',[0.526 1 0.526],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize);
       %overplot the objective part in another color
        h3 = parallelcoords(mPlot(:,1:(nO+1)),'color',[1 0.526 1],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize);
        
        %shouldn't be a need to overplot the optimum
        if 0
        for k=1:nOpts
            if (nargin>12) && (iOpt(k)>0)
                lLineWidth=2;
                %highlight optimal solution in another color
                h4 = parallelcoords(mPlot(iOpt(k),1:n),'color',[0 0 0],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize);
                h3 = parallelcoords(mPlot(iOpt(k),1:(nO+1)),'color',[0 0 0],'LineStyle','-','linewidth',lLineWidth,'Standardize',UseStandarize);
            end
        end
        end
    end
    
    %% Cleanup the plot
    
    set(plot2,'yLim',[ymin ymax],'Xticklabel',[],'Xtick',[],'fontsize',fontsize-4);

    %print the title + axis labels
    mytitle = sprintf('Near optimal region');
    %title(mytitle,'FontSize',fontsize);
    ylabel(plot2, AxisLabelsNew{1},'fontsize',fontsize);
    
    %duplicate axes to plot decision variable values at right
    ax2 = axes('Position',get(plot2,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'YColor',mColors(1,1,2,:));
       
    [ymin ymax];
                
    set(ax2,'yLim',[ymin ymax],'xLim',[1 n+0.5],'Xticklabel',[],'Xtick',[],'fontsize',fontsize-4);
    RightTicks = str2num(get(ax2,'YTickLabel'));
    class(RightTicks);
    RightTicks;
    ylabel(ax2,AxisLabelsNew{2},'fontsize',fontsize);

    
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
    

%plot the evolving range for x1 to x2

%startpt = 1.6;
%ruse = r;
%utangent = zeros(1,envd+1);
%ltangent = zeros(1,envd+1);


%utangent(2) = startpt;
%ltangent(2) = startpt;
%utangent(1)= c(2)*(utangent(2)-c(3)).^3+c(5)*(utangent(2)-c(6)).^2+c(7)*utangent(2)-c(4);
%ltangent(1)=utangent(1);

%for i=3:envd-1    
%    utangent(i) = sqrt(ruse^2 - sum(utangent(2:i-1).^2));
%    ltangent(i) = -sqrt(ruse^2 - sum(ltangent(2:i-1).^2));
    
    %check objective function bound
    
%           if i==2
%                ubnd =  c(2)*(utangent(i)-c(3)).^3+c(5)*(utangent(1)-c(6)).^2+c(7)*utangent(1)-c(4);
%                utangent(i) = max([ubnd utangent(i)]);
%                ltangent(i) = max([ubnd ltangent(i)]);
%           end

%end

%utangent;
%ltangent;

%h = plot(utangent,'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
%h = plot(ltangent,'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);

   
    % Determine the ticks on the right axis scale
    TicksRight = get(ax2,'ytick');
    TickSize = TicksRight(2)-TicksRight(1);
    
    %Rebuild the variable argument list to pass along through the controls
    
%    if strcmpi(AxisScales,'custom')
        vararg_curr = {'hWind' hWindReturn 'tolerance' Tolerance 'SubSpaceError' SubSpaceError 'NumSamples' NumSamples 'vFixed' vFixed 'Precision' Precision 'vFixedVals' vFixedVals ...
                'vStep' vStep 'fontsize' fontsize 'NumTicks' NumTicks 'TickMag' TickMag 'mActCat' mActCat ...
                'vGroup' vGroup 'mGroupData' mGroupData 'GroupToHighlight' iOpt 'SubSpaceGroup' SubSpaceGroup 'mColors' mColors 'mHighlightColor' mHighlightColor ...
                'AMat' AMat 'Brhs' Brhs 'cFunc' cFunc 'PlotPosition' PlotPosition 'PanelWidth' PanelWidth 'vObjLabels' vObjLabels ...
                'vXLabels' vXLabels 'vXLabelsShort' vXLabelsShort 'yAxisLabels' yAxisLabels 'AxisScales' AxisScales 'BaseAxis' BaseAxis 'mLims' mLims ...
                'TickSize' TickSize 'sGamsFile' sGamsFile 'StartTab' StartTab, 'GenerateType' GenType 'GenerateMethod' GenMethod 'HideSliders' HideSliders};

         %Also set app variables for the functon inputs and computations
         setappdata(hWindReturn,'varargs',vararg_curr);
         setappdata(hWindReturn,'mConvert',mTransformToOrig);
         setappdata(hWindReturn,'mDecs',mDecisions);
         setappdata(hWindReturn,'mObjs',mObjs);
         

%% Add controls to the right of the plot
% There are 5 basic elements
% 3 elements are frames for the Near Optimal, Sliders, and Display tabs that have lots
%     of controls
% 1 element is a frame and checkbox at the top right that allows the user
%     to hide all the tabs and widen the plot to cover the full figure width
% The final elements comprise the siders and axis check boxes on the plot
%     itself         
         
   if ShowControls == 1

    %define a local fontsize for the controls
    fontsizecntls = min([fontsize-4 12]);

    % Add the frames and controls
    %Frames

    %Frame to turn off/on all the other frames and controls
    ControlToggleFrame = uipanel('Title','','FontSize',fontsizecntls,...
                 'BackgroundColor','white',...
                 'Position',[xLeft+xWidth+lMarginPanel yBottom+yHeight+0.015 xPanelWidthNorm .045]);

    %Frame to contain the tabs         
    ControlFrame = uipanel('Title','','FontSize',fontsizecntls,...
                   'Position',[xLeft+xWidth+lMarginPanel 0.025 xPanelWidthNorm yBottom+yHeight],'BackgroundColor',get(hWindReturn,'Color'),'BorderType','none'); %

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
    
    
    %Box and controls for Generate Solutions box
    GenSols = uipanel('Title','Generate New Alternatives','FontSize',fontsizecntls-2,...
                 'BackgroundColor','white','parent',hTabs(2),...
                 'units','pixel','Position',[7 vPosTab(4)-243 267 137]); %[0.02 0.65 1-0.04 .2]);
             
    %Box and controls for Filter Existing Alternatives
    FiltSols = uipanel('Title','Filter Existing Alternatives','FontSize',fontsizecntls-2,...
                 'BackgroundColor','white','parent',hTabs(2),...
                 'units','normalized','Position',[0.02 0.01 1-0.04 .65-0.02]);
             
     set(FiltSols,'Units','pixel');
     vFiltPos = get(FiltSols,'position');
     set(FiltSols,'Units','normalized');


    %NearOptFrame = uipanel('Title','Near Optimal Controls','FontSize',fontsizecntls,...
    %             'BackgroundColor','white','parent',ControlFrame,...
    %             'Position',[xLeft+xWidth+0.08 0.025 .18 yBottom+yHeight-0.05]);

    %SliderValsFrame = uipanel('Title','Set Slider (Fixed) Values','FontSize',fontsizecntls,...
    %             'BackgroundColor','white','parent',ControlFrame,...
    %             'Position',[xLeft+xWidth+0.08 0.025 .18 yBottom+yHeight-0.025]);

    %DisplayFrame = uipanel('Title','Display Controls','FontSize',fontsizecntls,...
    %             'BackgroundColor','white','parent',ControlFrame,...
    %             'Position',[xLeft+xWidth+0.08 0.025 .18 yBottom+yHeight-0.025]);

    %define these control arrays here so the sliders callback can see them
    
    
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
        cbChecks(i) = uicontrol('Parent',Figure2,'Style', 'checkbox', 'Units','normalized','Position',[xLeft+(i-1)*xWidth/(n-0.5)-0.005 yBottom-0.03 0.02 0.02],'fontsize', fontsizecntls,'Value',vFixed(i),...
                            'Callback',{@AxisChecked,cbChecks,i});%
    end
    
    sCount = 1;
    maxSliders = 15;
    
    %Check All Axes Checkbox
    if sum(vFixed)==n
        AllValue = 1;
    else
        AllValue = 0;
    end        
        
    %Check box for hide/show slider control box
    %cbShowSliderFrame = uicontrol('Parent',hTabs(2),'Style', 'checkbox', 'String', 'Show slider panel','Position',[10 275 150 20],'fontsize', fontsizecntls,'Callback',{@ShowSliderFrame,hTabs(2)},'Value',1);    
    
    sSetBeyond = sprintf('When sliders are visible, check to set a value\non an axis beyond the range shown by the slider');
    
    cbHideSliders = uicontrol('Parent',FiltSols,'Style', 'checkbox', 'Value',1,'String', 'Hide sliders','Position',[8 vFiltPos(4)-45 106 20],'fontsize', fontsizecntls,'value',HideSliders); %,'Callback',{@HideSliders,sSlider}); [20 lTopSlider 165 20]
    cbCheckAll = uicontrol('Parent',FiltSols,'Style', 'checkbox', 'String', 'Check all axes','Position',[123 vFiltPos(4)-45 123 20],'fontsize', fontsizecntls,'Value',AllValue); %,'Callback',{@CheckAllBoxes,cbChecks},);% [20 lTopSlider-20 165 20]
    cbAllowSets = uicontrol('Parent',FiltSols,'Style', 'checkbox','ToolTipString',sSetBeyond, 'String', 'Set beyond extents','Position',[48 vFiltPos(4)-70 165 20],'fontsize', fontsizecntls,'Value',AllValue); %,'Callback',{@CheckAllBoxes,cbChecks},);% [20 lTopSlider-40 165 20]

    
    %Create the sliders
    
    %size(vMinPlot)
    %size(vMaxPlot)
    %size(vFixedValsPlotUnits)
    %size(vSliderStep)
    
    [sSlider,vShowSliders] = RenderSliders(vMinPlot,vMaxPlot,vFixedValsPlotUnits',vSliderStep);

    if 0
        if vMinResult(i) < vMaxResult(i)
            vShowSliders(i) = 1;
            
            vSliderStepValue = [vSliderStep(i) vSliderStep(i)];

            sSlider(i) = uicontrol('Style', 'slider',...
                'Min',vMinPlot(i),'Max',vMaxPlot(i),'Value',vFixedValsPlotUnits(i),'SliderStep',vSliderStepValue,'Visible','off',...
                'Units','normalized','Position', [xLeft+(i-1)*xWidth/(n-0.5)-0.005 yBottom+yHeight*(vMinPlot(i)-ymin)/(ymax-ymin) 0.01 yHeight*(vMaxPlot(i)-vMinPlot(i))/(ymax-ymin)]);
        else
            sSlider(i)=0;
        end    

    end
    
    for i=1:n
              
%        if vShowSliders(i) > 0 
            txtSliderValue(i) = uicontrol('Parent',FiltSols,'Style', 'edit', 'Position', [10 vFiltPos(4)-110-17*(sCount+1.75) 55 15],'fontsize', fontsizecntls-2, 'Callback', {@SetSliderValue,sSlider,mTransformToOrig,i});

            lblSliderValue(i) = uicontrol('Parent',FiltSols,'Style', 'text','String',vAxisLabelsAll{i}, 'Position', [70 vFiltPos(4)-110-17*(sCount+1.75) 110 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
            %set(txtSliderValue(i),'String',sprintf('%.2f',get(sSlider(i),'Value')));
            set(txtSliderValue(i),'String',sprintf('%.2f',ConvertPlot2DataUnits(get(sSlider(i),'Value'),mTransformToOrig(:,i))));
            sCount = sCount+1;
%        end
    end     
    %circle back and now set the callback function for the sliders
    for i=1:n
       set(sSlider(i),'Callback',{@SetTxtValues,txtSliderValue,mTransformToOrig,i})
       set(txtSliderValue(i),'Callback', {@SetSliderValue,sSlider,mTransformToOrig,i});
    end  

    %Controls to set all sliders to a particular solution
    % Set to specified record #
    lblSetAll = uicontrol('Parent',FiltSols,'Style', 'text','String','Set to:', 'Position', [10 vFiltPos(4)-95 44 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');%[10 lTopSlider-63 60 15]
    txtCurrRec = uicontrol('Parent',FiltSols,'Style', 'edit','Position', [55 vFiltPos(4)-95 45 15],'fontsize', fontsizecntls-2,...
              'Callback', {@UpdateSetToControls,hWindReturn,0,0});%[10 lTopSlider-63 60 15]  
    set(txtCurrRec,'String',lCurrRecord);
    
    sTotCount = sprintf('of %d',m);
    lblOfX = uicontrol('Parent',FiltSols,'Style', 'text','String',sTotCount, 'Position', [105 vFiltPos(4)-95 65 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');%[10 lTopSlider-63 60 15]

    %Buttons to advance records
    recFarLeft = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to the first record','String', '<<',...
            'Position', [55 vFiltPos(4)-115 20 15],'fontsize', fontsizecntls,....
            'Callback', {@UpdateSetToControls,hWindReturn,-1,2}); 
    recLeft = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to prior record','String', '<',...
            'Position', [75 vFiltPos(4)-115 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,-1,1});
    recRight = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'ToolTipString','... to next record','String', '>',...
            'Position', [95 vFiltPos(4)-115 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,1,1});
    recFarRight = uicontrol('Parent',FiltSols,'Style', 'pushbutton','ToolTipString','... to last record', 'String', '>>',...
            'Position', [115 vFiltPos(4)-115 20 15],'fontsize', fontsizecntls,...
            'Callback', {@UpdateSetToControls,hWindReturn,1,2});
    
    if lOptGroupCurr == 0
        sEnable = 'off';
    else
        sEnable = 'on';
    end
    
    GroupToHighlightButton = uicontrol('Parent',FiltSols,'Style', 'pushbutton', 'String', sprintf('Group to\nHighlight'),...
            'Position',[167 vFiltPos(4)-115 87 40] ,'fontsize', fontsizecntls,'enable',sEnable,...
            'Callback', {@SetSliderValsToRecord,lHighRecord});
            % 'Callback', {@SetAllSliderVals,sSlider,txtSliderValue,mTransformToOrig,2,mPlot(OptGroupInds,:)}); %[150 lTopSlider-65 70 20]

    %Tool tips
    sNearOpt = sprintf('Tolerance is the fraction of the optimal objective function value.\n Change from 1.0 to add a constraint to the underlying optimization problem\nto only consider solutions with objective functon values within\nthe specified tolerance.');
    sSubSpace = sprintf('Error is the fraction of a graph tick mark and used to define\n when a solution is within an error bound of the fixed solution\ndefined by checked axes and set values');
    
    %Label and Textbox for user to enter near-optimal tolerance
    NOToleranceLabel = uicontrol('Parent',hTabs(2),'Style', 'text','String','Near Optimal Tolerance (fraction of optimal):', 'Position', [10 lNearTopInteract-20 195 40],'fontsize', fontsizecntls,'BackgroundColor','white');
    NearOptTolerance = uicontrol('Parent',hTabs(2),'Style', 'edit','ToolTipString',sNearOpt, 'Position', [215 lNearTopInteract-20+10 45 20],'fontsize', fontsizecntls);
    set(NearOptTolerance,'String',Tolerance);
    
    lblSubSpaceError = uicontrol('Parent',hTabs(2),'Style', 'text','String','Subspace Error (fraction of tick):', 'Position', [10 lNearTopInteract-40 195 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    txtSubSpaceError = uicontrol('Parent',hTabs(2),'Style', 'edit','ToolTipString',sSubSpace, 'Position', [215 lNearTopInteract-40 45 20],'fontsize', fontsizecntls-2,'Callback',@TriggerSubSpaceUpdate);
    set(txtSubSpaceError,'String',SubSpaceError);
    

    
    %Drop-down box for Generate Type
    lTopGen = 90;
    GenerateType = 1;
    GenerateTypeLabel = uicontrol('Parent',GenSols,'Style', 'text','String','Type:','HorizontalAlignment','left', 'Position', [10 lTopGen 40 20],'fontsize', fontsizecntls-2,'BackgroundColor','white');
    cbGenerateType = uicontrol('Parent',GenSols,'Style','popup','String','One solution|Maximum extents|Random sample|Enumerate all (MIPs)','Position',[50 lTopGen 125 25],'fontsize',fontsizecntls-2,'Value',GenerateType,'Callback',{@SetFromGenBoxes});
    set(cbGenerateType,'value',GenType);
    
    GenerateMethodLabel = uicontrol('Parent',GenSols,'Style', 'text','String','Using:','HorizontalAlignment','left','Position', [10 lTopGen-30 40 20],'fontsize', fontsizecntls-2,'BackgroundColor','white');   
       
    sUsingTip = sprintf('Data - Query data already on the plot\nMatlab - Run Matlab linprog function using constraint system defined by\n          the parameters AMat, Brhs, and (possibly) cFunc\nGAMS - Run the specified GAMS file'); 
    cbGenerateMethod = uicontrol('Parent',GenSols,'Style','popup','ToolTipString',sUsingTip,'String','Data|MATLAB (LP matrix)|GAMS','Position',[50 lTopGen-30 125 25],'fontsize',fontsizecntls-2,'Value',GenMethod,'Callback',{@SetFromGenBoxes});
 
    
    lblNumSamples = uicontrol('Parent',GenSols,'Style', 'text','String','# Samples:', 'Position', [10 lTopGen-55 75 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    txtNumSamples = uicontrol('Parent',GenSols,'Style', 'edit', 'Position', [90 lTopGen-55 45 20],'fontsize', fontsizecntls-2);
    set(txtNumSamples,'String',num2str(NumSamples));
    %ResampleButton = uicontrol('Parent',hTabs(2),'Style', 'pushbutton', 'String', 'Resample',...
    %        'Position', [140 lTopSlider-30 90 20],'fontsize', fontsizecntls,'visible','off');
 
    lblGamsFile = uicontrol('Parent',GenSols,'Style', 'text','String','GAMS File:', 'Position', [10 lTopGen-80 75 15],'HorizontalAlignment','Left','fontsize', fontsizecntls-2,'BackgroundColor','white');
    sGamsFileTip = 'Full path and file name to GAMS file';
    txtGamsFile = uicontrol('Parent',GenSols,'Style', 'edit', 'ToolTipString',sGamsFileTip,'Position', [90 lTopGen-80 170 20],'fontsize', fontsizecntls-2);
    set(txtGamsFile,'String',sGamsFile);
    %RunGamsButton = uicontrol('Parent',hTabs(2),'Style', 'pushbutton', 'String', 'Run GAMS',...
    %        'Position', [180 lTopSlider-60 90 20],'fontsize', fontsizecntls,'visible','off');

    sGenSolsTip = sprintf('Generate new alternatives from the specifiednear-optimal tolerance\nand existing filtered alternatives');   
    GenerateButton = uicontrol('Parent',GenSols,'Style', 'pushbutton', 'String', 'Generate','ToolTipString',sGenSolsTip,...
            'Position', [185 lTopGen-25 75 40],'fontsize', fontsizecntls);
 
        
    
    %Check box for show/hide all controls (print view) -- goes above main frame
    cbShowControls = uicontrol('Parent',ControlToggleFrame,'Style','checkbox', 'String', 'Show all controls','Position',[5 7 145 20],'fontsize', fontsizecntls,'Callback',{@ShowAllControls,ControlFrame,[plot2 ax2],cbChecks,sSlider,hTabContainer,xPanelWidthNorm},'Value',1);

    %Controls for the File Tab
    
    % Save current data
    SaveButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Save Data',...
            'Position', [5 lTopLine 115 20],'fontsize', fontsizecntls,'callback',{@SaveFig,hWindReturn});

    % Buttons for plot/clear subspace
    PlotSubSpace = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Plot Subspace',...
            'Position', [5 lTopLine-50 115 20],'fontsize', fontsizecntls);
    ClearSubSpace = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Clear Subspace',...
            'Position', [120 lTopLine-50 120 20],'fontsize', fontsizecntls);        % Pushbutton string callback

    %Button to load pareto optimal solutions
    ParetoCheck = uicontrol('Parent',hTabs(1),'Style', 'checkbox', 'String', 'Add','Position',[100 lTopLine-88 60 20],'fontsize', fontsizecntls);
    ParetoButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load Pareto',...
            'Position', [10 lTopLine-88 90 20],'fontsize', fontsizecntls,...
            'Callback', {@LoadPareto,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep,sSlider,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,ParetoCheck});        % Pushbutton string callback

    %Button to permuate subregion
    TestButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Permute Subregion',...
            'Position', [28 lTopLine-110 128 20],'fontsize', fontsizecntls,...
            'Callback', {@TestAllRegion,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep,sSlider,mTransformToOrig,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,vGroup,1});
    LoadTestResults = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load Results',...
            'Position', [28 lTopLine-135 128 20],'fontsize', fontsizecntls,...
            'Callback', {@TestAllRegion,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep, sSlider,mTransformToOrig,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,vGroup,0});     
        
        
    %Button to plot the subregion
    PlotButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Plot Outer',...
            'Position', [28 lTopLine-160 128 20],'fontsize', fontsizecntls,...
            'Callback', {@RunButton,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,1,vararg_curr{:}});        % Pushbutton string callback
                                       % that calls a MATLAB function                                        
                                       
    %Button to reload gams output
    LoadGAMSButton = uicontrol('Parent',hTabs(1),'Style', 'pushbutton', 'String', 'Load GAMS Output',...
            'Position', [28 lTopLine-185 128 20],'fontsize', fontsizecntls,...
            'Callback', {@RunButton,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,0,vararg_curr{:}});        % Pushbutton string callback
        %'Callback', {@RunButton,hWindReturn,mResults,nO,NearOptTolerance,cbChecks,vFixedVals,vStep, sSlider,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,vGroup,0});        % Pushbutton string callback
                                       % that calls a MATLAB function
   
                                       
    %Controls for the Display Tab    

    %Controls for ramping color on traces along one axis
    CurrColorRamp = 0;
    CurrDirection = 1;
    sRampTip = sprintf('Ramp line colors on first checked axis\nin specified direction over specified # of classes\n(e.g., light to dark)');
    sDirectionTip = sprintf('Direction of color ramp:\nAscend: light to dark (small to large values)\nDescend: dark to light (small to large values)');
    cbRampColor = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', 'Ramp color','ToolTipString',sRampTip,'Position',[15 lTopColor 250 20],'fontsize', fontsizecntls,'Callback',{@RampColor},'Value',CurrColorRamp);
    lblNumClasses = uicontrol('Parent',hTabs(3),'Style', 'text', 'String','# Color classes:','Position',[10 lTopColor-1*25 130 20],'fontsize', fontsizecntls-2,'BackgroundColor','white','HorizontalAlignment','Left');
    txtNumClasses = uicontrol('Parent',hTabs(3),'Style', 'edit', 'Position',[110 lTopColor-1*25 30 20],'fontsize', fontsizecntls-2,'String','10');
    lblDirection = uicontrol('Parent',hTabs(3),'Style','text','String','Direction:','Position',[140 lTopColor-1*25 60 20],'fontsize', fontsizecntls-2,'BackgroundColor','White');
    cbDirection = uicontrol('Parent',hTabs(3),'Style','popup','ToolTipString',sDirectionTip,'String','Ascend|Descend','Position',[200 lTopColor-1*22 75 20],'fontsize', fontsizecntls-2,'Callback',{@RampColor},'Value',CurrDirection);
    
    lNumRowsForColorChecks = floor((lTopColor-30 -(lNearTop-3*20))/20);
    
    %Check box for show decision variable labels
    sShowGrouping = sprintf('Show a table below the axes labels\ncomprised of inputs from mGroupData\n that separates & colors decision axes\nby major and minor categories');
    sHideChecks = sprintf('Hide the check boxes below each axis');
    cbGroupDecisionLabels = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'ToolTipString',sShowGrouping, 'String', 'Show grouping labels','Position',[15 lNearTop-3*20 180 20],'fontsize', fontsizecntls,'Callback',{@GroupDecisionLabels,hTexts,hTabContainer,nO,plot2},'Value',0);
     %lblCheckAll = uicontrol('Parent',hTabs(1),'Style', 'text','String','Check all axes', 'HorizontalAlignment','Left','Position', [40 110 110 20],'fontsize', fontsizecntls,'BackgroundColor','white');
    cbHideChecks = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'ToolTipString',sHideChecks, 'String', 'Hide checkboxes','Position',[15 lNearTop-4*20 180 20],'fontsize', fontsizecntls,'Callback',{@HideCheckboxes,cbChecks,hTexts,hTabContainer,ymin-6*(ymax-ymin)/100,ymin-2*(ymax-ymin)/100,0.025,'normalized'});



    %Label and Textbox for move first checked axis left or right
    MoveAxisLabel = uicontrol('Parent',hTabs(3),'Style', 'text','String','Move checked axes:', 'HorizontalAlignment','center', 'Position', [20 lNearTop-6.75*20 100 40],'fontsize', fontsizecntls,'BackgroundColor','white');

    sFarLeft = ['... to the far left'];
    sFarRight = ['... to the far right'];
    sLeft = ['... one position to the left'];
    sRight = ['... one position to the right'];
    
    FarLeftButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton', 'ToolTipString',sFarLeft,'String', '<<',...
            'Position', [125 lNearTop-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults, nO, 1, -1,NearOptTolerance,txtSubSpaceError,cbChecks, sSlider,mTransformToOrig,vararg_out{:}});
        %     'Callback', {@MoveAxis,hWindReturn,mResults, nO, 1, -1,NearOptTolerance, cbChecks, vFixedVals, vStep, sSlider,vObjLabels,vXLabels, vXLabelsShort, yAxisLabels, fontsize, iOpt, mActCat, vGroup});
 
    LeftButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton', 'ToolTipString',sLeft,'String', '<',...
            'Position', [145 lNearTop-6.75*20+15 20 20],'fontsize', fontsizecntls);
        %    'Callback', {@MoveAxis,mResults,nO, 0, -1,NearOptTolerance,txtSubSpaceError, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});
    RightButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton', 'ToolTipString',sRight,'String', '>',...
            'Position', [165 lNearTop-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults,nO, 0, 1,NearOptTolerance,txtSubSpaceError, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});
    FarRightButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','ToolTipString',sFarRight, 'String', '>>',...
            'Position', [185 lNearTop-6.75*20+15 20 20],'fontsize', fontsizecntls);
         %   'Callback', {@MoveAxis,mResults,nO, 1, 1,NearOptTolerance,txtSubSpaceError, cbChecks,sSlider,mTransformToOrig,vararg_out{:}});

         
   %checkbox that determines how axes are re-ordered
    sReorder = sprintf('Reorder axes on the plot\nby the selected feature.');
    sByCat = '... by category specified in the input mActCat';
    sByFin = sprintf('... by values on the axis\n(e.g., first plot axes with singular (positive) values,\nnext axes with positive ranges,\n last axes that all are zero)');
    cbByCat = uicontrol('Parent',hTabs(3),'Style', 'checkbox','ToolTipString',sByCat, 'String', 'By Category','Position',[90 lNearTop-8.5*20+5 115 20],'fontsize', fontsizecntls);
    cbReorder = uicontrol('Parent',hTabs(3),'Style', 'checkbox','ToolTipString',sByFin, 'String', 'Defin. actions 1st','Position',[90 lNearTop-9.5*20+5 115 20],'fontsize', fontsizecntls);%    


    %Button to reorder axes
    ReorderButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','ToolTipString',sReorder,...
            'Position', [15 lNearTop-9.5*20 70 50],'fontsize', fontsizecntls,...
            'Callback', {@Reorder,mResults,nO,hWindReturn});
                    %'Callback', {@Reorder,mResults, nO,NearOptTolerance, cbChecks, sSlider,mTransformToOrig,cbReorder,cbByCat,vararg_out{:}});        % Pushbutton string callback

    %Wrap the text in ReorderButton
    cText = {'Reorder axes'};
    [outstring,newpos] = textwrap(ReorderButton,cText);
    set(ReorderButton,'String',outstring); %,'Position',newpos);        
                 

        
    %Button to reorder decision variables
    PruneButton = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','ToolTipString','Remove selected axes', ...  %'HorizontalAlignment','center',
            'Position', [15 lNearTop-12*20 70 50],'fontsize', fontsizecntls); %set the callback function after define SubSpaceError
        
    %Wrap the text in PruneButton
    cText = {'Remove axes'};
    [outstring,newpos] = textwrap(PruneButton,cText);
    set(PruneButton,'String',outstring); %,'Position',newpos);    
        
    %checkboxes and buttons to prune axes
    cbPruneChecked = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', 'Checked','Position',[90 lNearTop-11*20+5 115 20],'fontsize', fontsizecntls);
    cbPruneZeros = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', 'All zero','Position',[90 lNearTop-12*20+5 115 20],'fontsize', fontsizecntls);%    

    

    %Trace Group checkboxes to toggle their visibility
    if nU>0
        cbGroupChecks = zeros(nU,1);
        txtGroupOrders = zeros(nU,1);
        txtGroupNames = zeros(nU,1);
        for i=1:nU
            if i==lOptGroupCurr
                cColor = mHighlightColor;
            else
                lIndToUse = 1+mod(i-1,size(mColors,1));
                cColor = mColors(lIndToUse,1,3,:);
            end

            cbGroupChecks(i) =  uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', '','Position',[15 lNearTop-(16+i)*20+5 20 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'value',vShowGroup(i));
            txtGroupOrders(i) =  uicontrol('Parent',hTabs(3),'Style', 'edit','ToolTipString','Number to specify group order in plot','String', num2str(i),'Position',[35 lNearTop-(16+i)*20+5 28 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'BackgroundColor','white');
            txtGroupNames(i) =  uicontrol('Parent',hTabs(3),'Style', 'edit','ToolTipString','Label for group','String', vUnique(i),'Position',[65 lNearTop-(16+i)*20+5 145 20],'fontsize', fontsizecntls,'ForegroundColor',cColor,'BackgroundColor','white','HorizontalAlignment','left', ...
                                    'Callback',{@RenameGroup,i});
        end

        cbGroupsOneColor = uicontrol('Parent',hTabs(3),'Style', 'checkbox','ToolTipString','Plot all groups with a single color','String', 'All Groups Same Color','Position',[15 lNearTop-(17+nU)*20+5 195 20],'fontsize', fontsizecntls,'ForegroundColor',[0 0 0],'value',0);
        cbNoHighlight = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', 'No Highlight','Position',[15 lNearTop-(18+nU)*20+5 195 20],'fontsize', fontsizecntls,'ForegroundColor',[0 0 0],'value',0);


        ShowGroupsBut = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','String','Show Checked Groups', 'HorizontalAlignment','center', 'Position', [15 lNearTop-14*20+5 195 20],'fontsize', fontsizecntls, ...
              'Callback', {@ShowGroups,hWindReturn,[2 1]});
        RemoveGroupsBut = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','String','Remove Checked Groups', 'HorizontalAlignment','center', 'Position', [15 lNearTop-15*20+5 195 20],'fontsize', fontsizecntls, ...
              'Callback', {@RemoveGroups,mResults,nO,hWindReturn,2});
        sReorderGrp = sprintf('Reorder groups by number');
        ReorderGroupsBut = uicontrol('Parent',hTabs(3),'Style', 'pushbutton','ToolTipString',sReorderGrp,'String','Reorder Groups', 'HorizontalAlignment','center', 'Position', [15 lNearTop-16*20+5 195 20],'fontsize', fontsizecntls);
    else
        cbGroupChecks=[];
        txtGroupOrders=[];
        txtGroupNames=[];
    end
 
    %Build a cell array of handles to the controls we'll need to query in
    %callback functions. Organize by single-value text Box inputs, Axis
    %inputs, Group inputs, and other UI settings.
    control_handles = {'Tolerance' NearOptTolerance 'SubSpaceError' txtSubSpaceError 'NumSamples' txtNumSamples ...
                        'Tabs' hTabs 'AxisChecked' cbChecks 'Sliders' sSlider 'txtSliderValue' txtSliderValue ...
                        'GroupChecks' cbGroupChecks 'GroupOrders' txtGroupOrders 'GroupNames' txtGroupNames ...
                        'HideSliders' cbHideSliders 'ShowGroupDecisionLables' cbGroupDecisionLabels 'HideChecks' cbHideChecks 'ByCat' cbByCat 'Reorder' cbReorder 'PruneChecked' cbPruneChecked ...
                        'PruneZeros' cbPruneZeros 'AllGroupsSameColor' cbGroupsOneColor 'HighlightAsNormal' cbNoHighlight 'GamsFile' txtGamsFile ...
                        'ColorRamp' cbRampColor 'RampDirection' cbDirection 'NumClasses' txtNumClasses 'cbGenType' cbGenerateType 'cbGenMethod' cbGenerateMethod ...
                        'CurrentRecord' txtCurrRec 'TotalRecords' lblOfX 'GroupToHightlightBut' GroupToHighlightButton}; 

    %Also set as an app variable
    setappdata(hWindReturn,'hControls',control_handles);                
                    
    set(FarLeftButton,'Callback', {@RearrangeAxes,1,-1,hWindReturn,mResults,nO});
    set(LeftButton,'Callback', {@RearrangeAxes,0,-1,hWindReturn,mResults,nO});
    set(RightButton,'Callback', {@RearrangeAxes,0,1,hWindReturn,mResults,nO});
    set(FarRightButton,'Callback', {@RearrangeAxes,1,1,hWindReturn,mResults,nO});

    
    set(PlotSubSpace,'Callback', {@SubSpace,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,txtSubSpaceError,cbGroupChecks,1,vararg_curr{:}});        % Pushbutton string callback
    set(ClearSubSpace,'Callback', {@SubSpace,mResults,nO,NearOptTolerance,cbChecks,sSlider,mTransformToOrig,txtSubSpaceError,cbGroupChecks,0,vararg_curr{:}});        % Pushbutton string callback
    set(PruneButton,'Callback', {@RearrangeAxes,0,0,hWindReturn,mResults,nO});
    %set(ResampleButton,'Callback',{@Resample,NearOptTolerance,txtSubSpaceError,txtNumSamples,cbGroupChecks,cbChecks,sSlider,mTransformToOrig,mResults,nO,hWindReturn});
    %set(ResampleButton,'Callback',{@Resample,hWindReturn,control_handles,mTransformToOrig,mResults,nO});
    %set(RunGamsButton,'Callback',{@EnumerateWithGams,hWindReturn,control_handles,mTransformToOrig,mResults,nO});
    set(GenerateButton,'Callback',{@GenerateNewSols,hWindReturn,mTransformToOrig,mResults,nO});
    
    set(RemoveGroupsBut,'Callback',{@RemoveGroups,mResults,nO,hWindReturn,2});
    set(ReorderGroupsBut,'Callback',{@ReorderGroups,mResults,nO,hWindReturn})
    %New callback to allow calling nearoptplot2 at the end
    set(LoadTestResults,'Callback',{@TestAllRegion,hWindReturn,mResults,nO,control_handles,mTransformToOrig,0});
    
    set(cbHideSliders,'Callback',{@callHideSliders,sSlider});
    set(cbCheckAll,'Callback',{@CheckAllBoxes,cbChecks});%


    %set(PruneButton,'Callback', {@Prune,NearOptTolerance,txtSubSpaceError,cbChecks,sSlider,mTransformToOrig,cbPruneChecked,cbPruneZeros,mResults,nO,vararg_out{:}});        % Pushbutton string callback
    %Old version        'Callback', {@Prune,hWindReturn,mResults, nO,NearOptTolerance, cbChecks, vFixedVals, vStep, sSlider,vObjLabels,vXLabels, vXLabelsShort, yAxisLabels, fontsize, iOpt, mActCat,vGroup, cbPruneChecked, cbPruneZeros});        % Pushbutton string callback
    
    
else
    cbGroupChecks=[];
end

%Set Figure units to characters to allow multi-platform plotting
set(hWindReturn,'Units','characters');
%Set all other object units to normalized to allow multi-platform plotting
set(cbChecks,'Units','normalized');
set(hTexts,'Units','normalized');
for i=1:3
    set(get(hTabs(i),'Children'),'Units','normalized');
end
set(get(GenSols,'Children'),'Units','normalized');
set(txtSliderValue,'Units','normalized');
set(lblSliderValue,'Units','normalized');
set(cbGroupChecks,'Units','normalized');
set(txtGroupOrders,'Units','normalized');
set(txtGroupNames,'Units','normalized');

%% Callback functions for mouse events
% These events highlight lines on the plot

%Set the mouse call backs for moving the mouse and clicking to fix a value
set(hWindReturn,'WindowButtonMotionFcn', @MouseHoverCallback);
set(hWindReturn,'WindowButtonDownFcn', @MouseButtonDownCallback);
set(hWindReturn,'ResizeFcn',@ResizeCallback);


%Create a dummy over- parallel coordinate plot to highlight mouse-overed
%solutions

% Find out which groups are checked
gDexes = GetAllChecked(cbGroupChecks);

hold on
MouseOverPCPDec=parallelcoords([],'color',mColors(gDexes(1),1,3,:)); % GreenToMagentaRamp(1));
SubSpacePCPDec=parallelcoords([],'color',mColors(gDexes(1),1,2,:)); %GreenToMagentaRamp(3));

%Create the ramps 
[ColorRampPCPs,hColorChecks] = CreateEmptyColorRampPCPs(txtNumClasses);

if nO>0
    MouseOverPCPObj=parallelcoords([],'color',mColors(gDexes(1),2,3,:)); %GreenToMagentaRamp(16));
    SubSpacePCPObj=parallelcoords([],'color',mColors(gDexes(1),2,2,:)); %GreenToMagentaRamp(13));
end

textHdl = text('Color', 'black', 'VerticalAlign', 'Bottom');

%% Adding a Pareto inset plot
% 
% If there are pareto points and there are two objective, plot an inset x-y graph of them in standard
% cartesian coordinates in the upper left corner

if blHasPareto && (nO==2)
   hold off
   %first calculate where the inset figure should go and draw a box there.
   main_fig = findobj(Figure2,'Type','axes');
   ax=cell2mat(get(main_fig,'Position'));
   inset_size = 0.2;
   InnerPos = [.7*ax(1,1)+ax(1,3)-0.8*inset_size .9*ax(1,2)+ax(1,4)-inset_size 0.8*inset_size inset_size];
  
   FigureInset = figure('Color',[0.8 0.8 0.8]); 
   
   plot1 = subplot(1,1,1,'FontSize',fontsize-4,'Parent',FigureInset,'box','on','Color',[0.8 0.8 0.8]);
   hold on
   plot(plot1, mGroupSave(:,1),mGroupSave(:,2),'marker','x','MarkerSize',fontsize-4,'color',mColorsObjsGroup(2,:),'linestyle','none')
   xlabel(sLabel{1},'FontSize',fontsize-2,'color',[.737 0 0.737])
   ylabel(sLabel{2},'FontSize',fontsize-2,'color',[.737 0 0.737])
   %whitebg(FigureInset,[0.8 0.8 0.8]);
   %Inset the figure
   %inset(Figure2,FigureInset,0.25);
   
   inset_fig = findobj(FigureInset,'Type','axes');
   h_inset = copyobj(inset_fig,Figure2);
  
   set(h_inset,'Position', InnerPos);
  
   %outer_pos = get(h_inset,'OuterPosition');
   %rect = rectangle('Position',InnerPos,'EdgeColor',[0 0 0],'FaceColor',[0.8 0.8 0.8]);
   
   close(FigureInset);   
end

%% The Mouse Call back functions

    function ResizeCallback(src,evt)
        %Callback when the figure is resized. Grab the current data and
        %control settings and replot the figure (to update the positioning)
       
        ReorderGroups(0,0,mResults,nO,hWindReturn);
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
    
        axesHdl = plot2;
        %delete the existing mouse-overed parallelcoords plot and create a
        %new dummy one
        delete(MouseOverPCPDec);
        if nO>0
            delete(MouseOverPCPObj);
        end
        hold on
        
        %Find the groups to work with
        OrderedGroups = vUnique; %unique(vGroup);
        [gDexes] = GetAllChecked(cbGroupChecks);
        vLineStyleToUse = vLineStyle{gDexes(1)};
        
        %Use that group to see the colors
        MouseOverPCPDec = parallelcoords([],'color',mColors(gDexes(1),1,3,:)); %GreenToMagentaRamp(1));
        if nO>0
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
        % Read in the SubSpaceErrorVal
        if ShowControls==1
            [SubSpaceErrorVal, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'MouseMove',TickSize);

            %Find the rows of mResults that are checked and fixed
            [cRows] = ReturnRows(vGroup,OrderedGroups(gDexes),SubSpaceErrorValPU,mResults,cbChecks,vFixedVals,sSlider,mTransformToOrig);
            %ErrorTolPU = ConvertPlot2DataUnits(ErrorTolDU,mTransformToOrig(:,cI),1);

            %Search along the cI-th column of mResults to find values that within
            %the ErrorTol of the y-mouse position
            RowsToHighLight=cRows(abs(mPlot(cRows,cI)-mouseY) <= SubSpaceErrorValPU);
        else
            RowsToHighLight=[];
        end
        
        xrange = range(get(axesHdl, 'Xlim'));
        yrange = (get(axesHdl, 'Ylim'));
        
        if 1
            %set(textHdl, 'String', {['cI=', num2str(cI)], ['y=', num2str(mouseY)],['cGroupValue=', num2str(cGroupValue)], ['vGroup=', num2str(vGroup')],['cRows=', num2str(size(cRows))], ['R To High=', num2str(size(RowsToHighLight))]});
            %set(textHdl, 'String',{['SSEV=',num2str(SubSpaceErrorVal)],['SSEVPU=',num2str(SubSpaceErrorValPU)],['DV=',num2str(ConvertPlot2DataUnits(mouseY,mTransformToOrig(:,cI),0))],['PV=',num2str(mouseY)]});
            %set(textHdl, 'String',{[ThousandSep(ConvertPlot2DataUnits(mouseY,mTransformToOrig(:,cI),0))],[num2str(gDexes(1))]})
            set(textHdl, 'String',ThousandSep(ConvertPlot2DataUnits(mouseY,mTransformToOrig(:,cI),0)));
        else
            set(textHdl,'String','');
        end  

       % set(textHdl, 'String', {['cI=', num2str(cI)], ['y=', num2str(mouseY)],['cGroupValue=', num2str(cGroupValue)], ['cRows=', num2str(cRows)], ['R To High=', num2str(RowsToHighLight)]});
        set(textHdl, 'Position', [cI + 0.01*xrange, mouseY + 0.01*yrange])
        
        if ~isempty(RowsToHighLight)
            MouseOverPCPDec = parallelcoords(mPlot(RowsToHighLight,:),'color',mColors(gDexes(1),1,3,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize);
            hold on
            if nO>0
                MouseOverPCPObj = parallelcoords(mPlot(RowsToHighLight,1:nO+1),'color',mColors(gDexes(1),2,3,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize);
            end
           
        end
   end

    function UpdateHighlightedTraces
        % Updates the highlighted traces that are overplotted in darker
        % colors
        % 
        % The Source of traces depends on whether sliders are hidden or
        % visible:
        %   - Hidden: The current record (if specified)
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
        SubSpacePCPDec = parallelcoords([],'color',mColors(gDexes(1),1,2,:)); %GreenToMagentaRamp(3));
        if nO>0
            SubSpacePCPObj = parallelcoords([],'color',mColors(gDexes(1),2,2,:)); %GreenToMagentaRamp(13));
        end

        %Read whether the hide sliders check is selected
        %Query the relevant controls and structures
        [vParams,hControls] = ReadControls('UpdateHighlightTraces',hWindReturn);
        [sGamsFile,vFixed,vFixedVals,ToleranceValue,lCurrRecord,LimitSource,lHideSliderValue] = aGetField(vParams,{'sGamsFile' 'vFixed' 'vFixedVals' 'Tolerance' 'CurrentRecord' 'GenerateMethod' 'HideSliders'});
       
        %Error check validity of limit settings. If not valid, default to a
        %LimitSource setting of 1 (from Data records)     
        if (LimitSource==2) && (isempty(AMat) || isempty(Brhs))
            warning('Linear system not defined (AMat and/or Brhs). Continuing with plot data');
            set(cbGenerateMethod,'Value',1);
            LimitSource = 1;
        end
        if (LimitSource==3) && ~exist(sGamsFile,'file')
            warning('Gams file does not exist. Continuing with plot data');
            LimitSource = 1;
        end
       
        n = length(cbChecks);
        mResultsPlot = [];
 
        if (lHideSliderValue==1)
            %Sliders hidden, just highlight current record
            mResultsPlot = mPlot(lCurrRecord,:);
            
        elseif (LimitSource==1)
            %Sliders visible, draw from exists of existing data
            [mResCompact,mResultsData,cRows] = ReturnDataExtents(hWindReturn,mResults,nO,1);           
            mResultsPlot = mResultsData;
        elseif (LimitSource==2)       
            %Sliders visible, query min and max values using linear constraints defined in
            %AMat and brhs. Use maximum extent functions
            [AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,SubSpaceErrorValue] = UpdateLPMatrix(hWindReturn,nO);
            [mResultsValUse,mResCompact] = maxextentind(AMatNew,BrhsNew);
            for i=1:n
                mResultsPlot(:,i) = ConvertPlot2DataUnits(mResultsValUse(:,i),mTransformToOrig(:,i),1);
            end           
        else
            %Use GAMS to return the limits
            %Change the directory/folder to the one the GAMS File is in
            [sPath,sFileName,sFileExt] = fileparts(sGamsFile);            
            cd(sPath);

            RunMode = 2;
            ColInds = [1:n];      

            %Call GAMS to identify the maximum extents            
            [vObjs, mResultsInt, mResultsVal, uelsOut, vReturnFlag, mGamsStats, NumSolvs] = EnumNEIntSolsGams4(sGamsFile,vFixed(nO+1:nO+nD),vFixedVals(nO+1:nO+nD),vXLabelsShort,ToleranceValue,2,RunMode);

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
            mResCompact = [[min(vObjs(vReturnFlag==RunMode,:)); max(vObjs(vReturnFlag==RunMode,:))] [min(mResultsValUse);max(mResultsValUse)]];
            
            if RunMode==2
                ColInds = ColInds(vFixed==0);
            end
        end  
        
        %Show darkened lines for the current identified solutions (sub-set)
         if ~isempty(mResultsPlot)
            SubSpacePCPDec = parallelcoords(mResultsPlot,'color',mColors(gDexes(1),1,2,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize);
            hold on
            if nO>0
                SubSpacePCPObj = parallelcoords(mResultsPlot(:,1:nO+1),'color',mColors(gDexes(1),2,2,:),'LineStyle',vLineStyleToUse,'linewidth',1,'Standardize',UseStandarize);
            end
            
            NewLimitsPU = zeros(2,n);
            vCurrValsPU = zeros(1,n);

            if lHideSliderValue ~= 1
                for i=1:n
                    for j=1:2
                        NewLimitsPU(j,i) = ConvertPlot2DataUnits(mResCompact(j,i),mTransformToOrig(:,i),1);
                    end

                    vCurrValsPU(i) = ConvertPlot2DataUnits(vFixedVals(i),mTransformToOrig(:,i),1);
                end

                [sSliderOut,vShowSlider] = RenderSliders(NewLimitsPU(1,:),NewLimitsPU(2,:),vCurrValsPU,vSliderStep,sSliderUse);           
            end
         else
             SubSpacePCPDec = parallelcoords([],'color',mColors(gDexes(1),1,2,:));
             if nO>0
                SubSpacePCPObj = parallelcoords([],'color',mColors(gDexes(1),2,2,:));
             end
         end        
    end

    function SetSliderValue(hind,event,sSlider,mConvertFactors,i)
        %called when a txtSliderValue value is entered
        %If cbAllowSets is unchecked, checks entry is allowed (within slider range) and then sets the corresponding sSlider value
        %i is the index number of the txtSliderValue and sSlider
        %
        %Remember, txtValues are in original Data units and sSlider values are
        %in Plot units so we use mCovertFactors to convert

        lBeyondExtents = get(cbAllowSets,'Value');
        
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
        
        lHideSliderValue = get(cbHideSliders,'value');
        
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

    function CheckAllBoxes(cbCheck,event,cbTargetChecks)
        %called when user checks the cbCheckAll box
        %sets all the axes checkboxes cbChecks to the same value

        CheckValue = get(cbCheck,'Value');
        [m n] = size(cbTargetChecks);
        for i=1:n
            set(cbTargetChecks(i),'Value',CheckValue);
        end
        UpdateHighlightedTraces
    end

    function [vSlidersOut,vShowSlider] = RenderSliders(vMins,vMaxes,vFixedValsPU,vSliderStep,vSliders,Indexes,TotalInds) 
        % Renders the height of the sliders according to the specified
        % range in vMins (minimums) and vMaxes(maximumns). Hides sliders that have the same min and max
        %
        % INPUTS
        % vFixedValsPU = is the current setting of the slider in Plot Units
        % vSliderStep = is the step for the slider to move then clicking an
        % end (in plot units)
        % Optional vSliders is an array of object ids to the sliders. If
        % passed, the function will update. If omited, the function will
        % create new arrays
        % Optional Indexes are the original indexes of the sliders (if a
        % subset are passed) to position correctly
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
            error('RenderSliders: Paremeters are incompatible sizes, vMins: %d %d, vMaxes: %d %d, vVixedValsPU: %d %d, vSliderStep: %d %d',vMins,vMaxes,vFixedValsPU,vSliderStep)
            return
        end
        
        %Initialize variables
        vShowSlider = vMins < vMaxes;
        
        lHideSliderValue = get(cbHideSliders,'value');
        
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
            
            if lHideSliderValue==1
                %The checkbox value overrules everything else.
                sVis = 'off';
            end
            
            %Check the fixed value is within the limits. If not, set to the
            %nearest limit
            SliderSetValue = vFixedValsPU(i);
             
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
                    'Units','normalized','Position', [xLeft+(i+Indexes(i)-1)*xWidth/(TotalInds-0.5)-0.005 yBottom+yHeight*(vMins(i)-ymin)/(ymax-ymin) 0.01 yHeight*(vMaxes(i)+Delta-vMins(i))/(ymax-ymin)]);
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

    function TriggerSubSpaceUpdate(hind,event)
        UpdateHighlightedTraces
    end

    function RenameGroup(hind,event,i)
        % Takes the newly renamed i'th element in txtGroupNames text box
        % and passes it onto vGroup and vUnique
        NewName = get(txtGroupNames(i),'string');
        
        if iscell(vGroup)
            GrpRows = strcmpi(vGroup,vUnique(i));
        else
            GrpRows = vGroup==vUnique(i);
        end
        lNumRows = sum(GrpRows);
        
        vGroup(GrpRows) = repmat(NewName,lNumRows,1);
        %vGroup;
        vararg_curr = aSetField(vararg_curr,'vGroup',vGroup);
        setappdata(hWindReturn,'varargs',vararg_curr);
        %aGetField(vararg_out,'vGroup')
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
       
       hPCPs = cell(lNumClasses);
       cbColors = zeros(lNumClasses);
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

    function RampColor(hind,event)
       % If the check box is checked, ramps the colors of traces in the specified direction based on
       % values on the first checked axis by turning on a new layer
       % Otherwise, turns off the ramping layer
       
       % Direction => Ascend: lightest to darkest (smallest value to
       %                largest value; darkest layer on top/last)
       %              Descend: reverse -- darkest to lightest (largest
       %              value to smallest value; darkest layer on top/last)
       
       % Read the checkbox setting
       lChecked = get(cbRampColor,'value');
       lDirection = get(cbDirection,'value');
       
       % Delete the existing color ramp layer
       [ColorRampPCPs, hColorChecks] = DeleteColorRampPCPs(ColorRampPCPs,hColorChecks,txtNumClasses);
       
       hold on
       
       if lChecked == 0
           %turn off the color ramp layer and turn back on the original
           %layers
           ShowGroups(0,0,hWindReturn,[2 1]);
       else
           %turn off the prior layers and turn on the color ramp layer for the first checked axis
           
           %turn off the prior layers
           for l = 1:size(hPCGroup,1) 
               for j=1:2
                   for k=1:max(size(hPCGroup{l,j}))
                       set(hPCGroup{l,j}(k),'Visible','off');
                   end
               end
           end
           
           %Find the Axis and number of color classes
           [sFirst,lFirst] = GetFirst(cbChecks, {vObjLabels{:} vXLabels{:}}, 1);

           %Read in the number of color classes
           lNumClasses = str2num(get(txtNumClasses,'string'));

           if isempty(lNumClasses)
               warning('Invalid input for # Color Classes (%s). Defaulting to 10',get(txtNumClasses,'string'));
               lNumClasses=10;
           end
           
           ranges = [min(mPlot(:,lFirst)) max(mPlot(:,lFirst))];
               
           if ranges(1) == ranges(2)
               warning('RampColor','Min and Max on %s (%d-th axis) are the same. Can not ramp color. Select a different axis.',sFirst{:},lFirst)
               %turn off the color ramp layer      
               return
           end

           lNumClasses = round(lNumClasses);   

           %Build the color ramp with the specified number of classes
           
           %BaseRamp = OSUColorRamps('LightToDarkBlue7Step'); %reshape(mColors(1,1,:,1:end),3,3)
           BaseRamp = (OSUColorRamps('GreenToMagenta16Step'));
                      
           mRamp = ExpandColorRamp(lNumClasses,flipud(BaseRamp(1:7,:)),0,0);
           
           %Build the position matrix
           
           vOffset = mod([1:lNumClasses]-1,lNumRowsForColorChecks-1)+1;
           hOffset = floor(([1:lNumClasses]-1)/(lNumRowsForColorChecks-1));
                      
           if lDirection == 2
               cRange = [lNumClasses 1];
           else
               cRange = [1 lNumClasses];
               vOffset = fliplr(vOffset);
               hOffset = fliplr(hOffset);

           end

           %Reasign the data to the classes. Use linear interpolation along the
           %axis
           vGroupNew = round(interp1(ranges,cRange,mPlot(:,lFirst)));
           vGroupRanges = interp1(cRange,ranges,[1:(lNumClasses-1)/lNumClasses:lNumClasses]);
           
           %Cycle through the classes, create one plot and checkbox per color/group
           lWidth = 100;
           for i=1:lNumClasses
                ColorRampPCPs{i} = parallelcoords(mPlot(vGroupNew==i,:),'color',mRamp(i,:));
                
                vPosition = [25+(lWidth+5)*hOffset(i) lTopColor-30-(vOffset(i))*20 lWidth 20];
                %if (lNumClasses <= lNumRowsForColorChecks) || (i <= ceil(lNumClasses/2));
                %    %First left column
                %    vPosition = [25 lTopColor-30-(vOffset(i))*20 lWidth 20];
                %else
                    %second right column
                %    vPosition = [30+lWidth lTopColor-30-(i-ceil(lNumClasses/2))*20 lWidth 20];
                %end
                
                hColorChecks(i) = uicontrol('Parent',hTabs(3),'Style', 'checkbox', 'String', ['Class ',num2str(i)],'Position',vPosition, ...
                    'fontsize', fontsizecntls-2,'Callback',{@ToggleColor,i},'Value',1,'visible','on', ...
                    'string',sprintf('%.1f to %.1f',vGroupRanges(i),vGroupRanges(i+1)),'ForegroundColor',mRamp(i,:));
           end
       end       
    end

    function ShowGroups(hObj,event,hWindCurr,vLineThicknesses)
        % Shows the traces for the groups that are checked; hides traces for
        %       groups that are not checked

        % hWindCurr = Handle to the current Figure
        % vLineThicknesses = 2 element vector of line thickness (first for
        %      GroupToHighlight, second regular

        [vParams,hControls] = ReadControls('ShowGroups',hWindCurr);
        [mGroupInfo] = aGetField(vParams,{'mGroupData'});
        [cbGroupsOneColor, cbNoHighlight] = aGetField(hControls,{'AllGroupsSameColor' 'HighlightAsNormal'});
        %read in the check box values
        
        cAllSameColor = get(cbGroupsOneColor,'Value');
        cNoHighlight = get(cbNoHighlight,'Value');
        vGroupChecked = cell2mat(mGroupInfo(:,2));
        nU = length(vGroupChecked);
        
        for i=1:nU
            if vGroupChecked(i)==1
                strVis = 'on';
            else
                strVis = 'off';
            end

            if cAllSameColor==1
                vGroupToUse = 1;
            else
                vGroupToUse = i;
            end   

            mColorToUse = squeeze(mColors(vGroupToUse,1:2,1,:));
            lLineThick = vLineThicknesses(2);

            if (cNoHighlight==0) && (i==lOptGroupCurr)
                for j=1:2
                    mColorToUse(j,:) = mHighlightColor;
                end
                lLineThick = vLineThicknesses(1);            
            end

            %make the group visible
            for j=1:2
               for k=1:max(size(hPCGroup{i,j}))
                   set(hPCGroup{i,j}(k),'Visible',strVis,'Color',mColorToUse(j,:),'linewidth',lLineThick);
               end
            end
        end
        
        UpdateSetToControls(0,0,hWindCurr,0,0);        
    end

    function UpdateSetToControls(hind,but,hWindCurr,direction,magnitude)
        % Updates the SetTo controls (record number and number of records)
        % based on the currently selected groups and the direction and
        % magnitude to advance records
        %
        % direction: -1 (left), +1 (right), 0 (none)
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

if 0
    get(plot2,'yLim')
    get(plot2,'ytick')
    get(plot2,'yticklabel')
    
    get(ax2,'yLim')
    get(ax2,'ytick')
    get(ax2,'yticklabel')
    mResults
    mPlot
end

hold off
end

%% Callback Functions for Controls on Plot

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
                set(hChildren(i),'Visible',Visibility);
            case 'uipanel'
               SetVisibilityAll(hChildren(i),Visibility);
        end
    end
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
    
    %xLeft = 0.100;
    %xWidth = 0.625;
    %yBottom = 0.3068;
    %yHeight = 0.6182;
    %xWidthAddWOControls = 0.18; %extra width to add without controls
    
    nAxes = max(size(hAxes));
    nControls = max(size(cChecks));
    nSlides = max(size(sSlides));
    %nLines = max(size(hLines));
    
    xWidthOld = xWidth;
    cVal = get(hind,'value');
    if cVal==0 %make all the frames invisible
        strVis = 'off';         
        xWidth=xWidth + xWidthAddWOControls; 
        
    else
        strVis = 'on';
        xWidth = xWidth-xWidthAddWOControls;
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
        i;
        vCheckPos(1)=xLeft+(i-1)*xWidth/(nControls-0.5)-0.005;
        
        set(cChecks(i),'Position',vCheckPos);
        if sSlides(i)~=0
            vSliderPos = get(sSlides(i),'Position');
            vSliderPos(1)=xLeft+(i-1)*xWidth/(nSlides-0.5)-0.005;
            set(sSlides(i),'Position',vSliderPos);
        end
    end
     
    xPosTable = get(hTabContainer,'Position');
    xPosTable(3) = xWidth*(nControls-0.5)/(nControls-1+0.5);
    set(hTabContainer,'Position',xPosTable);
    
    
    %for i=1:nLines
    %    if hLines(i)>0
    %        xPos = get(hLines(i),'x');

    %        [i xPos];
    %        xPos = xLeft+(xPos-xLeft)*xWidth/xWidthOld;
    %        xPos;
    %        set(hLines(i),'x',xPos);
    %    end
    %end
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
        %[num2str(i),'. ', fields{i}]
        %strcmpi(inArray,fields{i})
        %circshift(strcmpi(inArray,fields{i})',1)
        %inArray{circshift(strcmpi(inArray,fields{i})',1)}
        
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
end

function GroupDecisionLabels(hind,event,hTexts,hTabContainer, nO, hAxis)
    % toggles showing the decision variable axes labels in group mode with the
    % table of hTabContainer (checked)
    %
    % hTexts are the handles to the text labels in list mode
    % hTabContainer is the handle to the java element showing the table of
    % grouped decision variables

    % nO is an offset for the first nO elements of hTexts to leave on
    % regardless of the toggle (i.e., labels for the objective functions) 
    
    % hAxis is handle to the axis to which the Text Labels refer
    
    % also recalculates the position for the hTabContainer table just below
    % the hTexts
    
    cVal = get(hind,'value');
    cText= max(size(hTexts));
    %cGroup=max(size(hTextGroup));
    %cLines=max(size(hLines));
    
    if cVal==0 %make the hTexts visible and groups invisible
        sText = 'on';
        gText = 'off';
    else %keep the hTexts visible and make the groups visible
        sText = 'on';
        gText = 'on';
    end
    
    %mTextExtent = zeros(cText,4);
    mTextPos = zeros(cText,3);
    
    for i=nO+1:cText
        set(hTexts(i),'Visible',sText);
        %strCurrUnits = get(hTexts(i),'Units');
       % set(hTexts(i),'Units','normalized');
       % mTextExtent(i,:) = get(hTexts(i),'Extent');
        mTextPos(i,:) = get(hTexts(i),'Position');
        %set(hTexts(i),'Units',strCurrUnits);
    end
    
    yBottomTable = 0.025;
    aPosition = get(hAxis,'Position');
    yHeight = aPosition(4);
    yBottom = aPosition(2);
    
    %mMaxExtent = min(mTextExtent(nO+1:cText,2));
    mMaxPos = min(mTextPos(nO+1:cText,2));
    
    set(hTabContainer,'Visible',gText);
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

function [SubSpaceErrorVal, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,strCallFunction,TickSize)
    %returns the SubSpace Error Value in the checkbox. Checks that the
    %value is numerical and returns the absolute value
    %If the optional parameter TickSize is passed, additionally calculates
    %and returns SubSpaceErrorVal in Plot Units.
    
    SubSpaceErrorVal = get(txtSubSpaceError,'String');

    %error checking on SubSpaceErrorValue
    if ~isnumeric(str2num((SubSpaceErrorVal)))
            warning(['nearoptmo2: ',strCallFunction],'SubSpaceError is not a numerical value. Defaulting to zero to continue')
            SubSpaceErrorVal = 0;
    else
            SubSpaceErrorVal= abs(str2num(SubSpaceErrorVal));
    end
    
    if nargin>2
        SubSpaceErrorValPU=SubSpaceErrorVal*TickSize;
    else
        SubSpaceErrorValPU = NaN;
    end
end

function [vSliderValuesDU vSliderSteps vSliderValuesPU] = ReadSliderValues(sSliders,vFixedValues,mConvert)
    %returns the vector of sSliders values. Where a sSlider is not defined,
    %instead substitutes the predefined fixed value.
    %Note vFixedValues are provided in Data units, whiles sliders are in
    %plot units.
    %Outputs slider values in both data units (vSliderValuesDU) and plot
    %units (vSliderValuesPU)
    %Use the coefficients mConvert to convert between the two.
    
    vSliderValuesDU = vFixedValues;
    vSliderValuesPU = vSliderValuesDU;
    vSliderValuesPU(:)=0;
    vSliderSteps = vFixedValues;
    n = max(size(vSliderValuesDU));
    for i=1:n
        if sSliders(i)>0
            vSliderValuesPU(i) = get(sSliders(i),'Value');
            sSlideSteps = get(sSliders(i),'SliderStep');
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

function [outargs,hControls] = ReadControls(CallFunction,hWindCurr)
    %Queries select entries on hControls and updates corresponding values
    %in varargin and returns as outray %varargou
    
    % Get the handles for the controls we will need
    hControls = getappdata(hWindCurr,'hControls');
    varargin = getappdata(hWindCurr,'varargs');
    mConvert = getappdata(hWindCurr,'mConvert');
    
    [txtTolerance, txtNumSamples, cbChecks, txtSubSpaceError, sSliders, txtGroupOrders,txtGroupNames,cbGroupChecks,txtGamsFile,cTabs,txtCurrRec,cbGenType,cbGenMethod,cbHideSliders] = aGetField(hControls,{'Tolerance','NumSamples','AxisChecked','SubSpaceError','Sliders','GroupOrders','GroupNames','GroupChecks','GamsFile','Tabs','CurrentRecord','cbGenType','cbGenMethod' 'HideSliders'});
    %Group data
    mGroupData = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
    %Text boxes
    ToleranceValue = get(txtTolerance,'String');
    sGamsFile = get(txtGamsFile,'String');   
    if ~strcmpi(sGamsFile,'') && (exist(sGamsFile,'file') == 0)
        %File does not exist
        error(sprintf('nearoptplotmo2 (%s) - File %s does not exist',CallFunction,sGamsFile));
    end
    
    NumSamples = GetCheckSubSpaceError(txtNumSamples,CallFunction);
    %Axes check boxes
    vFixed = ReadCheckboxValues(cbChecks);
    %Start Tabs
    StartTab = GetAllValues(cTabs,'visible','on');
    %Current Record
    lCurrRecord = str2num(get(txtCurrRec,'String'));
    %Generate Type and method
    GenType = get(cbGenType,'value');
    GenMethod = get(cbGenMethod,'value');
    %Hide sliders
    HideSliders = get(cbHideSliders,'value');
        
    %Read in the variable arguments
    [vFixedVals TickSize] = aGetField(varargin,{'vFixedVals' 'TickSize'});   
    
    [SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,CallFunction,TickSize);
    [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    %Update the values in the field list
    outargs = aSetField(varargin,'Tolerance',ToleranceValue,'SubSpaceError',SubSpaceErrorValue,'SubSpaceErrorPU',SubSpaceErrorValPU,'NumSamples',NumSamples,'vFixed',vFixed,'vFixedVals',vFixedVals,'mGroupData',mGroupData,'sGamsFile',sGamsFile,'mConvert',mConvert,'StartTab',StartTab,'CurrentRecord',lCurrRecord,'GenerateType',GenType,'GenerateMethod',GenMethod,'HideSliders',HideSliders);
    %varargout = varargout_temp;
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

function [gNewOrder] = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks)
    % Takes the input on txtGroupOrders and txtGroupNames boxers all the group info (txtGroupNames, cbGroupChecks)
    % Returns a cell matrix of the newly sorted Group Names and Checked
    % values.

    %Read in the new order
    lNumGroups = length(txtGroupOrders);
    gOrder = cell(lNumGroups,3);
    for i=1:lNumGroups
        %i
%        if strcmpi(class(get(txtGroupOrders(i),'string')),'char')
%            gOrder(i,1) = {get(txtGroupOrders(i),'string')};
%            gOrder(i,2) = {str2num(get(txtGroupNames(i),'string'))};         
%        else
            gOrder(i,1) = {get(txtGroupOrders(i),'string')};
            gOrder(i,2) = get(txtGroupNames(i),'string');
%        end
        gOrder{i,3} = get(cbGroupChecks(i),'value');
    end
    
    %gOrder;

    [gSort, gIC] = sortrows(gOrder,1);
    
    %gSort;

    %Reassign the sorted variables
    gNewOrder = gSort(:,[2:3]);
end

function RemoveGroups(hObj,event,mResults,nO,hWindCurr,lAction)
    %Reassigns checked groups based on the lAction selected
    % 1 - re-assigns all checked groups to the first checked group
    % 2 - deletes the checked groups (permanently removes)
    
    %varargin_start = getappdata(hWindCurr,'varargs');
    [varargin_fresh,hControls] = ReadControls('RemoveGroups',hWindCurr);
    [vGroup iOpt mGroupData] =  aGetField(varargin_fresh,{'vGroup' 'GroupToHighlight' 'mGroupData'});

    [cbGroupChecks] = aGetField(hControls,{'GroupChecks'});
  
    vChecked = cell2mat(mGroupData(:,2));
    vUniques = mGroupData(:,1);

    %[iOpt] =  aGetField(varargin,{'GroupToHighlight'});

    %[vGroup iOpt] =  aGetField(,{'vGroup' 'GroupToHighlight'});
    m = size(vGroup,1);

    nG = length(vUniques);
    %Checked = ReadCheckboxValues(cbGroupChecks)

    vGroupNew = vGroup;
    
    switch lAction
        
        case 1 %re-assign all checked groups to the first checked group
    
            [gTemp, fInd] = GetFirst(cbGroupChecks,vUniques,0);
            gInds = [1:nG];
            gIndKeep = gInds((vChecked==0) | (gInds==fInd));
            
            gInds = [fInd+1:nG];  


            if iscell(vUniques)
                gRevertTo = vUniques{fInd};
            else
                gRevertTo = vUniques(fInd);
            end

            for i=gInds
                if vChecked(i)
                     if iscell(vGroupNew)
                         for j=1:m
                            if strcmpi(vGroupNew{j},vUniques(i))
                                vGroupNew{j}=gRevertTo;
                            end
                         end
                       else
                           vRevert= vGroupNew==vUniques(i);
                           vGroupNew(vRevert)=gRevertTo;
                     end
                end
            end

            mResultsNew = mResults;
            mGroupDataNew = mGroupData(gIndKeep,:);

        case 2 %delete the checked groups
            %Remove those rows from vGroup and mResults
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
            
    end 
    
    %Revert to default of show all retained groups
    for i=1:size(mGroupDataNew,1)
        mGroupDataNew{i,2} = 1;
    end
            
    varargout = aSetField(varargin_fresh,'vGroup',vGroupNew,'mGroupData',mGroupDataNew);
    [mObjs, mDecs]=SplitMatrix(mResultsNew,nO);
    nearoptplotmo2(mObjs, mDecs,varargout{:}); 
end

%function ReorderGroups(hObj,event,mResults,nO,txtGroupOrders,txtGroupNames,cbGroupChecks,varargin)
function ReorderGroups(hObj,event,mResults,nO,hWindCurr)
    %Reorders the groups by the entries in txtGroupOrders
    
    [vararg_out,hControls] = ReadControls('ReorderGroups',hWindCurr);
    %varargin_start = getappdata(hWindCurr,'varargs');
    %vararg_out = ReadControls('RemoveGroups',hControls,varargin_start{:});
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
    % Get the handles for the group controls we will need
    %[txtGroupOrders,txtGroupNames,cbGroupChecks] = aGetField(hControls,{'GroupOrders','GroupNames','GroupChecks'});
  
    %Update the groups
    %mGroupDataSort = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
    
    %Work on the grouping. The new group will appear with a child name just below the
    %first checked group in the resorted list
    
    %Find the first checked group
    [gValue, gDex] = GetFirst(cbGroupChecks,mGroupDataSort(:,1),1);

    if iscell(vGroup)          
        if nargin == 2
            % Add the specified name as a new full name below the first checked
            % group
            sNewGroupName = dNewGroup;
        else
           %auto generate the name from the 1st checked name. Append the
           %name part
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

    mGroupDataNew = [mGroupDataSort(:,1) mGroupDataSort(:,2) gSortStart(:); {sNewGroupName} {1} {gSortVal}];
    mGroupDataNewSort = sortrows(mGroupDataNew,3);
    
end
                                    
%function Reorder(hObj,event,mMatrix, nObjs, txtTolerance, cbChecks, sSliders,mConvert,cbReorder, cbByCat, varargin) %#ok<INUSL>
function Reorder(hObj,event,mMatrix, nObjs,hWindCurr)
  %  txtTolerance, cbChecks, sSliders,mConvert,cbReorder, cbByCat, varargin)
% Called to reorder the decision variable columns of mMatrix according to the dynamic ranges of the decision variables
    % If cbReorder is checked then put the constant non-zero decision variables first, dyanmic range variables second, and zero-valued variables third.
    % If cbReorder is NOT checked: put the variables with dynamic range
    % first, constant non-zero decision variables second, and zero-valued decision variables third
    
    % If cbByCat is checked then do the organization by the overarching
    % category type in vActCats (i.e., all 1's first, 2's second, etc...)
    
    % assumes first nObjs columns of mMatrix are objective function values
    % and are not part of the reordering
    
    % Also reorders the corresponding Fixed Axes in cbChecks, labels and
    % abreviations in vXLables and vXLabelsShort, and vActCats
    
    
    [max(mMatrix(:,2:end)); min(mMatrix(:,2:end))];
    
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
        
        %count columns in each category
        mActCat;
        cats=unique(mActCat(:,1))
        if strcmp(class(mActCat),'cell')
            mActCat = cell2mat(mActCat(:,1));
        end
        if strcmp(class(cats),'cell')
            cats = cell2mat(cats);
        end        
        CatCounts = histc(mActCat(:,1),cats);
        nCats = max(size(CatCounts));
       
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
    
    %nearoptplotmo(hWind, mObjs, mDecs, Tolerance, vNewFixed, vNewFixedVal,vNewStep, vObjLabels,vNewXLabels, vNewXLabelsShort,yAxisLabels, fontsize, iOpt, mNewActCat, vGroup);
end

%function RearrangeAxes(hObj,event,FarDir,Direction,txtTolerance,txtSubSpaceError,cbChecks,sSliders,mConvert,cbChecked,cbZero,mMatrix,nObjs,varargin) %#ok<INUSL>
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
    % cbChecked - checked (value of 1) and Direction of 0 == Prune/take out all the checked axes (columns in the matrix mMatrix and associated data).
    % cbZeros is checked and Direction of  0 = Prune/take out all zero-valued axes
    
    % Repositioning is only within the objective function or decision
    % variable axes -- i.e. can't move an objective function axes into and
    % among decision variable axes
    
    % nObjs = number of objective functions in mMatrix (first nObjs columns of mMatrix)
    
    % Also reorders/prunes the corresponding Fixed Axes in cbChecks, vObjLabels, and labels and
    % abreviations in vXLables and vXLabelsShort and all other columnar
    % data to facillitate replotting the chart
    
    % Step 1. Figure out new order of columns (from Direction, FarDir, cbChecked, cbZeros)
    % Step 2. Reorder all the columnar data (mMatrix, vFixedVals, vStep, vXLabels, vXLabelsShort, mActCat, vObjLabels, mLims, AMat, Brhs, cFunc
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
    ObjInds = [1:nObjs];
    
    if Direction==0 %Prune checked/zeroed axes
        %first identify which columns to {Prune (==1) and retain (==0)
        %start by assuming we keep them all (==0)
        vPruneCheck = zeros(1,n);
        vPruneZero = zeros(1,n);

        if iChecked>0
            vPruneCheck = vFixed;
            sum(vPruneCheck(1:nObjs),2);
            %if (sum(vPruneCheck(1:nObjs),2) == nObjs) && (nObjs > 0)
            %    error('nearoptplotmo2:Prune','Must keep at least one objective function checked');
            %    return
            %end
        end

        if iZero>0
            vMax = max(mMatrix);
            vMin = min(mMatrix);

            vPruneZero = (vMax == 0) .* (vMin==0);
        end

        vPrune = max([vPruneCheck; vPruneZero]);  
        vKeep = AxInds(vPrune==0);
        vKeepDec=DecInds(vPrune(nObjs+1:n)==0);
        vKeepObj=ObjInds(vPrune(1:nObjs)==0);
        vRemove = AxInds(vPrune==1);
        
        lPrune = sum(vPrune);
        lPruneObj = sum(vPrune(1:nObjs));
        lNewObj = nObjs-lPruneObj;
        lNewDec = nD - sum(vPrune(nObjs+1:n));
        
    else
        %Move checked axes left or right
            %calcualte the destination column
        vCheckDec =  DecInds(vFixed(nObjs+1:n)==1);
        vCheckObj = ObjInds(vFixed(1:nObjs)==1);
        vNoCheckDec = DecInds(vFixed(nObjs+1:n)==0);
        vNoCheckObj = ObjInds(vFixed(1:nObjs)==0);
        
        if FarDir==1           
            %goto the extremes
            if Direction==-1 %left means checked axes come first
                vKeepDec = [vCheckDec vNoCheckDec];
                vKeepObj = [vCheckObj vNoCheckObj];
            else %right means checked axes come last
                vKeepDec = [vNoCheckDec vCheckDec];
                vKeepObj = [vNoCheckObj vCheckObj]; 
            end
        else
            %move only one away in direction     
            [vDummy, vKeepDec] = ShiftNonZero(vFixed(DecInds+nObjs),Direction);
            [vDummy, vKeepObj] = ShiftNonZero(vFixed(ObjInds),Direction);
            
            %add back in indexes for columns that were not selected
        end
        
        if nObjs>0
            vKeep = [vKeepObj vKeepDec+nObjs];
        else
            vKeep = vKeepDec+nObjs;
        end
        
        lNewDec = nD;
        lNewObj = nObjs;
    end
    
    %vKeep;
    %vKeepObj;
    %vKeepDec;
    
    %Read in the varargin parameters that will be maniputed
    [vFixedVal vStep vXLabels vXLabelsShort mActCat vObjLabels mLims AMat Brhs cFunc BaseAxis Tolerance SubSpaceError] = aGetField(varargin_rrax,{'vFixedVals' 'vStep','vXLabels','vXLabelsShort','mActCat','vObjLabels','mLims','AMat','Brhs','cFunc','BaseAxis','Tolerance','SubSpaceError'});
    
    %[vFixedVal vSliderSteps vFixedValPU] = ReadSliderValues(sSliders,vFixedVal,mConvert); 
    %SubSpaceError = GetCheckSubSpaceError(txtSubSpaceError,'Prune');
 
    [mS,nS] = size(mActCat);
    size(vFixedVal);
    size(vFixed);
    
    size(vObjLabels);
    size(vXLabels);
    
    
    vObjLabels;
    vXLabels;
    
    if nObjs==0
        vLabelsAll = vXLabels';
    else
        vLabelsAll = {vObjLabels{:} vXLabels{:}};
    end
    
    size(vLabelsAll)
        
    %retain the columns to keep
    [mMatrixN, vFixedN, vFixedValN, vStepN,vLabelsAllN, mLimsN] = ReorderCols(vKeep,mMatrix,vFixed,vFixedVal',vStep,vLabelsAll,mLims);
    
    %vKeepDec
    %vXLabelsShort
    %mActCat
    vXLabelsShortUse = vXLabelsShort;
    
    if size(vKeepDec) ~= size(vXLabelsShort)
        vXLabelsShortUse = vXLabelsShort';
    end
    
    [vXLabelsShortN,mActCatN] = ReorderCols(vKeepDec,vXLabelsShortUse,mActCat');

    mActCatN = mActCatN';
    
    if ~isempty(AMat)
        AMatN = ReorderCols(vKeepDec,AMat);
        
        %Update the right hand side constraints. Substract off from the right-hand side contraint valus the fixed values
        %for the columns to be removed
        [nD lNewDec];
    
        if lNewDec ~= nD
            BrhsN = Brhs - AMat*(vPrune(nObjs+1:n).*vFixedVal(nObjs+1:n)');
        else
            BrhsN = Brhs;
        end
    else
        AMatN=AMat;
        BrhsN = Brhs;
    end
    if ~isempty(cFunc)
        cFuncTemp = ReorderCols(vKeepDec,cFunc);
        %now work on the rows (object functions
        cFuncN = ReorderCols(vKeep(1:nObjs),cFuncTemp');
        cFuncN = cFuncN';
    else
        cFuncN = cFunc;
    end
    
    %Work on the BaseAxis
    BaseAxisN = BaseAxis;
    if isempty(find(vKeepObj==BaseAxis(1)))
        BaseAxisN(1)=min(lNewObj,1);
    else
        BaseAxisN(1) = find(vKeepObj==BaseAxis(1));
    end
    
    if isempty(find(vKeepDec==BaseAxis(2)))
        BaseAxisN(2)=1;
    else
        BaseAxisN(2) = find(vKeepDec==BaseAxis(2));
    end
    
    [mObjs, mDecs]=SplitMatrix(mMatrixN,lNewObj);
    [vObjLabelsN, vXLabelsN]=SplitMatrix(vLabelsAllN,lNewObj);
        
    %Set parameters that were manipulated to pass via varargin
    varargout = aSetField(varargin_rrax,'tolerance',Tolerance,'vFixed',vFixedN','vFixedVals',vFixedValN','vStep',vStepN,'vObjLabels',vObjLabelsN','vXLabels',vXLabelsN', ...
                'vXLabelsShort',vXLabelsShortN,'mActCat',mActCatN,'SubSpaceError',SubSpaceError,'BaseAxis',BaseAxisN,'mLims',mLimsN,'AMat',AMatN,'Brhs',BrhsN,'cFunc',cFuncN);
    
    nearoptplotmo2(mObjs, mDecs,varargout{:});
       
    %nearoptplotmo(hWind, mObjs, mDecs, Tolerance, vFixedN, vFixedValN,vStepN, vObjLabels,vXLabelsN, vXLabelsShortN,yAxisLabels, fontsize, iOpt, mActCatN, vGroup);
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

function [vGroupSubSet] = ReturnRows(vGroup,vGroupSelect,SubSpaceErrorValPU,mMatrix,cbChecks,vFixedVals,sSliders,mConvert,vFixed)
    % Finds and returns a vector of the row numbers in vGroup that have the
    % value vGroupSelect
    %
    % If optional parameters mMatrix, cbChecks, sSliders,
    % are passed, then only returns the rows in mMatrix whose columns have
    % values within SubSpaceError tolerance for the checked columns.
    %
    % Remember, mMatrix, SubSpaceError, and vFixedVals are in Data units and sSlider
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
    % txtSubSpaceError = handle to a textbox with the error tolerance
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
    
    %SubSpaceErrorVal = GetCheckSubSpaceError(txtSubSpaceError,'ReturnRows');    
    %SubSpaceErrorVal;
    %class(SubSpaceErrorVal);
      
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
        
        %Convert SubSpaceError to Data Units for each column
        SubSpaceErrorDU = zeros(1,n2);
        for i=1:n2
            SubSpaceErrorDU(i) = ConvertPlot2DataUnits(SubSpaceErrorValPU,mConvert(:,i));
        end
        vGroupSubSetWithError = [];
        
        % loop through the rows of mMatrix
        for i=vGroupSubSet'            
               resid = abs(mMatrix(i,:).*vFixed - vFixed.*vFixedVals');
               %[i resid SubSpaceErrorValPU]
               
               %resid
               %SubSpaceErrorValPU

               if (sum(resid <= SubSpaceErrorDU)==n2)
                   vGroupSubSetWithError = [vGroupSubSetWithError; i];
               end
        end
        vGroupSubSet = vGroupSubSetWithError;
    end 
end

function [AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,SubSpaceErrorValue] = UpdateLPMatrix(hWindCurr,nO)       
    % Use the checked axes, slider set values, and Subspace Error control setting to define a new linear system of equations representing the feasible area 
    % that is defined by the original linear program components (AMat, Brhs, cFunc) and the current fixed value settings.
    %
    % This new system of equations is returned at AMatNew, BrhsNew, and
    % cFuncNew
    %
    % They are build as follows:
    % For decision variable axess, this is simply adding two constraints to the AMat:
    %       1) Decision Variable Value <= Fixed Value + Error
    %       2) Decision Variable Value >= Fixed Value - Error
    % For objective functions, also add two constraints to AMat
    %       1) c(X) <= Fixed Value + Error
    %       2) c(X) >= Fixed Value - Error (where c(X) is pulled from the
    %                       coefficiencts in cFunc.
    %  
    %  As an example:
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
    %  like AMat, Brhs, cFunc
    %  hControls = cell array of handles to controls on the figure to
    %  pull out current settings like fixed (set) axis, set values,
    %  Error value, etc.
    %  mConvert = matrix of axis conversion factors from plot to data units
    %  nO = number of objectives (columns) in the data set
    %
    % OUTPUTS
    %  AMatNew = updated a-matrix defining the constraint coefficients
    %   (m+2*d x nF) where m=original number of rows/constraints,
    %   d=number of fixed deicision variables + fixed objective function
    %   values, and nF = number of decision variables
    %  BrhsNew = updated right hand side constraint coefficients (m+2*d
    %         x 1)
    %  cFunc = nO x nF maxtrix of objective function coefficients for
    %       the nO objectives
    %  cFuncFree = version of cFunc with columns removed for decision variable values that fixed (vFixed == 1 and Error=0)
    %  vFixed = 1 x n vector of binaries with a value of 1 indicating the axis is fixed
    %  vFixedVals = 1 x n vector of values specified the fixed value
    %  SubSpaceErrorValue = string value of the Sub Space Error control

    % Get the handles for the figure controls we will need
    [varargin_fresh,hControls] = ReadControls('UpdateLPMatrix',hWindCurr);
    [vFixed, vFixedVals, TickSize, AMat, Brhs, cFunc, SubSpaceErrorValue, SubSpaceErrorValPU, mConvert] = aGetField(varargin_fresh,{'vFixed' 'vFixedVals' 'TickSize' 'AMat' 'Brhs' 'cFunc' 'SubSpaceError' 'SubSpaceErrorPU' 'mConvert'});   

    %[SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'UpdateLPMatrix',TickSize);
    %[vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);

    %[cbChecks, txtSubSpaceError, sSliders] = aGetField(hControls,{'AxisChecked','SubSpaceError','Sliders'});

    %size(cbChecks)
    %vFixed = ReadCheckboxValues(cbChecks);
    [mF nF] = size(vFixed);
    n = nF;
    %Calc number of decision variables
    nD = nF-nO;
    dInds = [1:nD];
    oInds = [1:nO];

    %Read in the variable arguments
    %varargin_fresh = getappdata(hWindCurr,'varargs');
    %[vFixedVals TickSize AMat Brhs cFunc, SubSpaceErrorValPU] = aGetField(varargin_fresh,{'vFixedVals' 'TickSize' 'AMat' 'Brhs' 'cFunc' 'SubSpaceErrorPU'});   

    %[SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'UpdateLPMatrix',TickSize);
    %[vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);

    %Build the new constraint matrix that includes the fixed decision and
    %objective function values
    if isempty(AMat) || isempty(Brhs)
        error('UpdateLPMatrix: Can not resample. Need to provide AMat and/or Brhs parameters as inputs to represent the constraint system [AMat][X] <= Brhs of the underlying optimization problem.')
    end

    AMatNew = AMat;
    BrhsNew = Brhs;

    dToAdd = dInds(vFixed(nO+1:n) ~= 0);
    oToAdd = oInds(vFixed(1:nO) ~= 0);

    if SubSpaceErrorValPU==0
       %We need to remove the variable from the constraint set and use a
       %modified (reduced) matrix set
       blIsReduced = 1;

       BrhsNew = BrhsNew - AMat(:,vFixed(nO+1:end)==1)*vFixedVals(vFixed(nO+1:end)==1);
       AMatNew = AMat(:,vFixed(nO+1:end)==0);
       cFuncNew = cFunc(:,vFixed(nO+1:end)==0);
    else
        blIsReduced = 0;
        for i=dToAdd
            cRow = circshift(eye(nD,1),i-1)';
            ErrorValDU = ConvertPlot2DataUnits(SubSpaceErrorValPU,mConvert(:,nO+i));
            vFixedVals(i+nO);
            AMatNew = [AMatNew;cRow;-cRow];
            BrhsNew = [BrhsNew; vFixedVals(i+nO) + ErrorValDU; - vFixedVals(i+nO) + ErrorValDU];
        end
        cFuncNew = cFunc;
    end

    if (sum(vFixed(1:nO)) > 0) && (isempty(cFunc) || size(cFunc,1) < sum(vFixed(1:nO)))
        warning('Need to specify parameter cFunc to fix objective function values or not enough cFuncs specified. Continuing without')
    else
        for i=oToAdd
            ErrorValDU = ConvertPlot2DataUnits(SubSpaceErrorValPU,mConvert(:,i));
            AMatNew = [AMatNew;cFuncNew(i,:);-cFuncNew(i,:)];
            BrhsNew = [BrhsNew; vFixedVals(i) + ErrorValDU; - vFixedVals(i) + ErrorValDU];     
        end
    end

    %Build the matrix of free (not fixed) objective functions and their coefficients;
    %and pivot it
    if isempty(cFuncNew)
        cFuncFree = [];
    else
        cFuncFree = cFunc(vFixed(1:nO)==0,:)';    
    end

    %Output the matricies to check
    %size(cFunc);
    %size(cFuncFree);
    %cFuncFree;
    %size([AMat Brhs])
    %size([AMatNew BrhsNew])

    %BrhsNew(size(AMat,1)+1:end);
    %AMatNew(size(AMat,1)+1:end,:);
end

function [mResCompact,mResultsData,cRowsRet] = ReturnDataExtents(hCurr,mData,nO,blStrictExtent)
     %Query the max and min values from the data records in mData
     % defined by cbChecks (fixed axes) within SubSpaceError of
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
    [SubSpaceErrorVal,SubSpaceErrorValPU,mGroupDataSort,vGroup,vFixedVals,mConvert] = aGetField(varargin_fresh,{'SubSpaceError','SubSpaceErrorPU','mGroupData','vGroup','vFixedVals','mConvert'});
    
    %Read in select controls
    [cbChecks,sSlider] = aGetField(hControls,{'AxisChecked','Sliders'});
    
    [m,n] = size(mData);
    
    %[SubSpaceErrorVal, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'ReturnDataExtents',TickSize);
    ColInds = [1:n];

    %Find the groups to work with
    selGroups = mGroupDataSort(cell2mat(mGroupDataSort(:,2))==1,1); %unique(vGroup);
    cRows = ReturnRows(vGroup,selGroups,SubSpaceErrorValPU,mData,cbChecks,vFixedVals,sSlider,mConvert);
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
            cRowsCurr = ReturnRows(vGroup,selGroups,SubSpaceErrorValPU,mData,cbChecks,vFixedVals,sSlider,mConvert,vFixedUse);
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

%function Resample(hWind,event,txtTolerance,txtSubSpaceError,txtNumSamples,cbGroupChecks,cbChecks,sSliders,mConvert,mData,sControls,nO,hWindCurr)
function [mObjsNew, mDecsNew] = Resample(hWind,event,hWindCurr,mData,nO)
    %Re-sample an additional specified number of solutions from the sub-space defined by the checked axes and fixed values associated with those
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
    mConvert = aGetField(varargin_fresh,{'mConvert'});
    %Update the constraint matrix
    [AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,SubSpaceErrorValue] = UpdateLPMatrix(hWindCurr,nO);
    %Query the remaining controls
    %[txtTolerance txtNumSamples, txtGroupOrders,txtGroupNames,cbGroupChecks,txtGamsFile] = aGetField(hControls,{'Tolerance','NumSamples','GroupOrders','GroupNames','GroupChecks','GamsFile'});
    %ToleranceValue = get(txtTolerance,'String');
    %sGamsFile = get(txtGamsFile,'String');
    %NumSamples = GetCheckSubSpaceError(txtNumSamples,'Resample');   
    
    [ToleranceValue, sGamsFile, NumSamples, vGroup, SubSpaceGroup, mActCat, mGroupDataSort, AMat, Brhs, SubSpaceErrorValPU] = aGetField(varargin_fresh,{'tolerance' 'sGamsFile' 'NumSamples' 'vGroup' 'SubSpaceGroup' 'mActCat' 'mGroupData' 'AMat' 'Brhs' 'SubSpaceErrorPU'});
    
    [m n] = size(mData);
    nD = n-nO;
    %Read in the variable arguments
    %varargin_fresh = getappdata(hWindCurr,'varargs');
    %[vGroup SubSpaceGroup mActCat] = aGetField(varargin_fresh,{'vGroup' 'SubSpaceGroup' 'mActCat'});   
    
    %Update the groups
    %mGroupDataSort = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
        
   
    if 0
    
    [txtTolerance txtNumSamples, cbChecks, txtSubSpaceError, sSliders, txtGroupOrders,txtGroupNames,cbGroupChecks,txtGamsFile] = aGetField(hControls,{'Tolerance','NumSamples','AxisChecked','SubSpaceError','Sliders','GroupOrders','GroupNames','GroupChecks','GamsFile'});

    ToleranceValue = get(txtTolerance,'String');
    sGamsFile = get(txtGamsFile,'String');
    NumSamples = GetCheckSubSpaceError(txtNumSamples,'Resample');
    
    %size(cbChecks)
    vFixed = ReadCheckboxValues(cbChecks);
    [mF nF] = size(vFixed);
    [m n] = size(mData);
    
    %Calc number of decision variables
    nD = n-nO;
    dInds = [1:nD];
    oInds = [1:nO];
        
    %Read in the variable arguments
    varargin_fresh = getappdata(hWindCurr,'varargs');
    [vFixedVals vGroup SubSpaceGroup TickSize AMat Brhs cFunc mActCat] = aGetField(varargin_fresh,{'vFixedVals' 'vGroup' 'SubSpaceGroup' 'TickSize' 'AMat' 'Brhs' 'cFunc' 'mActCat'});   
    
    %Update the groups
    mGroupDataSort = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
    
    [SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'Resample',TickSize);
    [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    %Build the new constraint matrix that includes the fixed decision and
    %objective function values
    if isempty(AMat) || isempty(Brhs)
        error('Resample: Can not resample. Need to provide AMat and/or Bhrs parameters as inputs to represent the constraint system.')
    end
    
    AMatNew = AMat;
    BrhsNew = Brhs;
    
    dToAdd = dInds(vFixed(nO+1:n) ~= 0);
    oToAdd = oInds(vFixed(1:nO) ~= 0);
    
    if SubSpaceErrorValPU==0
       %We need to remove the variable from the constraint set and use a
       %modified (reduced) matrix set
       blIsReduced = 1;
       
       BrhsNew = BrhsNew - AMat(:,vFixed(nO+1:end)==1)*vFixedVals(vFixed(nO+1:end)==1);
       AMatNew = AMat(:,vFixed(nO+1:end)==0);
       cFuncNew = cFunc(:,vFixed(nO+1:end)==0);
    else
        blIsReduced = 0;
        for i=dToAdd
            cRow = circshift(eye(nD,1),i-1)';
            ErrorValDU = ConvertPlot2DataUnits(SubSpaceErrorValPU,mConvert(:,nO+i));
            vFixedVals(i+nO);
            AMatNew = [AMatNew;cRow;-cRow];
            BrhsNew = [BrhsNew; vFixedVals(i+nO) + ErrorValDU; - vFixedVals(i+nO) + ErrorValDU];
        end
        cFuncNew = cFunc;
    end
    
    if (sum(vFixed(1:nO)) > 0) && (isempty(cFunc) || size(cFunc,1) < sum(vFixed(1:nO)))
        warning('Need to specify parameter cFunc to fix objective function values or not enough cFuncs specified. Continuing without')
    else
        for i=oToAdd
            ErrorValDU = ConvertPlotToDataUnits(SubSpaceErrorValPU,mConvert(:,i));
            AMatNew = [AMatNew;cFuncNew(i,:);-cFuncNew(i,:)];
            BrhsNew = [BrhsNew; vFixedVals(i) + ErrorValDU; - vFixedVals(i) + ErrorValDU];     
        end
    end
    
    %Build the matrix of free (not fixed) objective functions and their coefficients;
    %and pivot it
    if isempty(cFuncNew)
        cFuncFree = [];
    else
        cFuncFree = cFunc(vFixed(1:nO)==0,:)';    
    end
        
    %Output the matricies to check
    %size(cFunc);
    %size(cFuncFree);
    %cFuncFree;
    %size([AMat Brhs])
    %size([AMatNew BrhsNew])
    
    %BrhsNew(size(AMat,1)+1:end);
    %AMatNew(size(AMat,1)+1:end,:);
    end
    
    %Check extents of matrixes
    [mExtOrig, mExtCompactOrig] = maxextentind(AMat,Brhs);
    [mExtNew,mExtCompactNew] = maxextentind(AMatNew,BrhsNew);
      
    %Resample
    NewSols = stratgibbs(NumSamples,AMatNew,BrhsNew,struct('matformat','reduce','lincombo',cFuncFree));
    
    [mNS nNS] = size(NewSols);
    
    %B. Work on the decision variables
    
    if sum(vFixed(nO+1:end)) > 0
        %Then must add fixed variable values back in for fixed variables
        NewSolsFull = zeros(mNS,nD);
        NewSolsFull(:,vFixed(nO+1:end)==1) = vFixedVals(vFixed(nO+1:end)==1);
        NewSolsFull(:,vFixed(nO+1:end)==0) = NewSols;
        mExtComNF = zeros(nD,2);
        mExtComNF(vFixed(nO+1:end)==1,:) = vFixedVals(vFixed(nO+1:end)==1);
        mExtComNF(vFixed(nO+1:end)==0,:) = mExtCompactNew;
    else
        % No need to change
        NewSolsFull = NewSols; 
        mExtComNF = mExtCompactNew;
    end
    
    mDecsNew = NewSolsFull;
    
    %blIsReduced;
    %size(mDecisions);
    %size(NewSolsFull);
    
    %mDecisionsNew = [mDecisions; NewSolsFull];
    
    %[size(AMatNew); size(BrhsNew); size(cFuncFree); size(AMat); size(Brhs); size(cFunc); size(NewSols); size(NewSolsFull)]

    [size([1:n-nO]');size(mExtCompactOrig); size(mExtComNF)]
    
    %Print out new decision ranges
    %sprintf('New Decision Ranges')
    ['# Orig Extents  New Extents']
    [[1:nD]' mExtCompactOrig mExtComNF]
    [num2cell([1:nD]') mActCat num2cell(mExtComNF)]
    
    %C. Work on the Objective functions; have to calculate objective function values for the new solutions
    if nO>0
        [size(NewSolsFull); size(cFunc')]
        mObjsNewSols = NewSolsFull*cFunc';
        sprintf('New Objective Function Range')
        %mObjsNew = [mObjs; mObjsNewSols];
        mObjsNew = mObjsNewSols;
        
        %Print out the range of new sampled obj function values
        [min(mObjsNewSols); max(mObjsNewSols)]
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
    [txtTolerance txtNumSamples, cbChecks, txtSubSpaceError, sSliders, txtGroupOrders,txtGroupNames,cbGroupChecks,txtGamsFile] = aGetField(hControls,{'Tolerance','NumSamples','AxisChecked','SubSpaceError','Sliders','GroupOrders','GroupNames','GroupChecks','GamsFile'});

    
    %ToleranceValue = get(txtTolerance,'String');
    %sGamsFile = get(txtGamsFile,'String');
    %NumSamples = GetCheckSubSpaceError(txtNumSamples,'Resample');
    
    %size(cbChecks)
    %vFixed = ReadCheckboxValues(cbChecks);
    [mF nF] = size(vFixed);
    [m n] = size(mData);
    
    %Calc number of decision variables
    nD = n-nO;
    dInds = [1:nD];
    oInds = [1:nO];
        
    %Read in the variable arguments
    %varargin_fresh = getappdata(hWindCurr,'varargs');
    [vFixedVals,vGroup,SubSpaceGroup,TickSize,vXLabelsShort,mGroupDataSort,SubSpaceErrorValue] = aGetField(varargin_fresh,{'vFixedVals' 'vGroup' 'SubSpaceGroup' 'TickSize' 'vXLabelsShort' 'mGroupData' 'SubSpaceError'});   
    
    %Update the groups
    %mGroupDataSort = SortGroupsByTextOrder(txtGroupOrders,txtGroupNames,cbGroupChecks);
    
    %[SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'Resample',TickSize);
    %[vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    %Change the directory/folder to the one the GAMS File is in
    if exist(sGamsFile,'file') == 0
        error(['nearoptplotmo2:EnumerateWithGAMS--the file ',sGamsFile,' does not exist.']);
        return
    end
    [sPath,sFileName,sFileExt] = fileparts(sGamsFile);
    cd(sPath);
    
    %[vXLabelsShort' num2cell(vFixed(nO+1:end)') num2cell(vFixedVals(nO+1:end))]
 
    %Pass the values and enumerate the integer solutions
    [mObjs, mResultsInt, mResultsVal, uelsOut, ReturnFlag, GamsStats, NumSolvs] = EnumNEIntSolsGams4(sGamsFile,vFixed(nO+1:end),vFixedVals(nO+1:end),vXLabelsShort,str2num(ToleranceValue),2,lRunMode);
    
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
    [sprintf('%d solves made\n%d new solutions added\n%d errors',NumSolvs,NumSolvs-numErrors,numErrors)]

    if numErrors>0
        ['Solutions with errors']
        [ReturnFlag(ReturnFlag<0) GamsStats(ReturnFlag<0,:) mResultsInt(ReturnFlag<0,:)]
    end        
end  


function GenerateNewSols(hWind,event,hWindCurr,mConvert,mData,nO)
    % Generates new solutions based on the selctions in the Generate
    % Solutions Box on the Interact Tab. Adds as a new group to the plot.

    %Read the controls
    [varargin_fresh,hControls] = ReadControls('GenerateNewSols',hWindCurr);    
    [vGroup,lGenUse,lGenType] = aGetField(varargin_fresh,{'vGroup' 'GenerateMethod' 'GenerateType'});   
    
    cTypeSuffix = {'Sing' 'Ext' 'Resample' 'Enum'};
    cUseSuffix = {'Dat' 'Mat' 'Gams'};
    
    cTypeSuffixFull = {'Single Solution' 'Maximum Extents' 'Resample' 'Eunumerate'};
    cUseSuffixFull = {'Data' 'Matlab' 'Gams'};
    
    if lGenUse==1
        %Data
        if lGenType~=2
            warning(['No function defined for GenType #', num2str(lGenType), '. Ignoring']);
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
        [AMatNew,BrhsNew,cFunc,cFuncFree,vFixed,vFixedVals,SubSpaceErrorValue] = UpdateLPMatrix(hWindCurr,nO);
       
        if lGenType==1
            %Single solution, optimize
            [mDecNew, mObjNew, exitflag, minoutput] = linprog(cFunc,AMatNew,BrhsNew,[],[],[],[],[],struct('maxiter',1000,'Display', 'off'));
            
        elseif lGenType==2
            %Maximum extents   
            [mDecNew,mResCompact] = maxextentind(AMatNew,BrhsNew);
            
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

        %A. Split the original data
        [mObjs,mDecisions]=SplitMatrix(mData,nO);

        %B. Work on the decision variables
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
        
        sprintf('%s from %s: %d solutions added',cTypeSuffixFull{lGenType},cUseSuffixFull{lGenUse},mNS)
        % Update the variables
        varargout = aSetField(varargin_fresh,'vGroup',[vGroup;vGroupAdd],'mGroupData',mGroupDataNewSort);  

        nearoptplotmo2(mObjsAll,mDecisionsAll,varargout{:});
    else
        h = errordlg('No new solutions generated','Solution Error');

        uiwait(h); 
    end
end

function SubSpace(hind,event,mMatrix,nObjs,txtTolerance,cbChecks,sSliders,mConvert,txtSubSpaceError,cbGroupChecks,iAction,varargin)
    %called when user presses the Plot/Hide sub space buttons
    %
    % iAction==1 : plot the subspace defined by the checked axes and
    % fixed/slider values for those axes. Does this by assigning a child identifier to the 1st checked vGroup
    % for each row in mMatrix that meets the conditions
    %
    % iAction==0 : remove the subspace. Does this by re-assigning each
    % vGroup entry with the child identifier back to a reverted identifier (1st checked vGroup)
    %
    
    %Collect the values from the controls
    ToleranceValue = get(txtTolerance,'String');
       
    vFixed = ReadCheckboxValues(cbChecks);
    [mF nF] = size(vFixed);
    %vFixedActions = vFixed((nObjs+1):n);
    
    %Read in the variable arguments
    [vFixedVals vGroup iOpt SubSpaceGroup TickSize] = aGetField(varargin,{'vFixedVals' 'vGroup' 'GroupToHighlight' 'SubSpaceGroup' 'TickSize'});   
    
    [SubSpaceErrorValue, SubSpaceErrorValPU] = GetCheckSubSpaceError(txtSubSpaceError,'SubSpace',TickSize);
    
    vGroup;
    SubSpaceGroup;
    
    vFixedVals;
    
    %fprintf('Before\n');
    %vFixedVals
       
    [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    [m n] = size(mMatrix);
    
    if size(vGroup,1) ~= m
        vGroup = ones(1,m);
    end
        
    %determine which vGroup is the orginal/base group
    [vUniques, iC, uIndex, iCBeg, iVBeg] = GroupOrderAppear(vGroup);
    
    vUniques;
    iCBeg;
    
    GroupRevert = vUniques(1);
    %find the first checked group
    gDexes = GetAllChecked(cbGroupChecks);
    GroupRevert = vUniques(gDexes(1));
    
    if 0
    if size(vUniques,1)>1
        for i=1:size(cbGroupChecks,1)
            if get(cbGroupChecks(i),'Value')==1
                GroupRevert = get(cbGroupChecks(i),'string');
                if isnumeric(vUniques)
                    GroupRevert = str2num(GroupRevert);
                end
                break
            end
        end
    end
    end
    GroupRevert;
    
    if iscell(vUniques)
        SubSpaceGroupNew = sprintf('%s_SubSpace%d',vUniques{gDexes(1)},size(vUniques,1)+1);
    else
        SubSpaceGroupNew = vUniques(gDexes(1))+(size(vUniques,1)+1)/10;
    end  
    
    if iAction==1
       %Find rows that are in the GroupRevert subspace defined by the checks and error
       %tolerance
       [vNewRows] = ReturnRows(vGroup,GroupRevert,SubSpaceErrorValPU,mMatrix,cbChecks,vFixedVals,sSliders,mConvert);
       size(vNewRows);
       repmat(SubSpaceGroupNew,size(vNewRows,1),1);
       if iscell(vGroup)
           
            vGroup{vNewRows} = repmat(SubSpaceGroupNew,size(vNewRows,1),1);
       else
            vGroup(vNewRows) = repmat(SubSpaceGroupNew,size(vNewRows,1),1);
       end
    else
       %Revert SubSpaceGroup rows to 
       [vNewRows] = ReturnRows(vGroup,SubSpaceGroup);
       if iscell(vGroup)
            vGroup{vNewRows} = GroupRevert;
       else
            vGroup(vNewRows) = GroupRevert;
       end
       SubSpaceGroupNew = 0;
    end
    
    
    vGroup;
    [mObjs mActVolumes] = SplitMatrix(mMatrix,nObjs);
    
    %split out group statistics
    if 0
    uGroup = unique(vGroup);
    uN = size(uGroup,1);
    
    gResults = zeros(uN,3);
    gResults(:,1) = uGroup;
    
    for i=1:uN
        vGMembers = find(vGroup == uGroup(i));
        gResults(i,2) = min(mMatrix(vGMembers,1));
        gResults(i,3) = max(mMatrix(vGMembers,1));
    end
        
    sprintf('Subspace objective function range\nGroup  Min     Max\n');
    gResults;
    end
    
    vGroup;
    SubSpaceGroupNew;
    %Update the parameters in the variable argument list;
    varargout = aSetField(varargin,'tolerance',ToleranceValue,'vFixed',vFixed,'vFixedVals',vFixedVals,'vGroup',vGroup,'SubSpaceError',SubSpaceErrorValue,'SubSpaceGroup',SubSpaceGroupNew);
    nearoptplotmo2(mObjs(:,1:nObjs), mActVolumes, varargout{:});
 end

function RunButton(hind,event,mMatrix,nObjs,txtTolerance,cbChecks,sSliders,mConvert,iRunGAMS,varargin)
%function RunButton(hind,event,hWind,mMatrix,nObjs,txtTolerance,cbChecks,vFixedVals,vStep,sSliders,mConvert,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,vGroup,iRunGAMS)
    %called when user presses the near optimal region button
    %Reads in Tolerance and check box values and writes to a text file GAMS
    %can read. Then calls GAMS (if iRunGAMS==1) and reads gams outputs back
    %into the figure
    
    %Collect the values from the controls
    ToleranceValue = get(txtTolerance,'String');
    vFixed = ReadCheckboxValues(cbChecks);
    [m n] = size(vFixed);
    vFixedActions = vFixed((nObjs+1):n);

    %Read in the variable arguments
    [vFixedVals vStep vXLabelsShort ] = aGetField(varargin,{'vFixedVals' 'vStep' 'vXLabelsShort'});   

    
    [vFixedVals vSliderVals vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    vMaxResult = max(mMatrix(:,(nObjs+1):n));
    vMinResult = min(mMatrix(:,(nObjs+1):n));
    
    mActCat;
    iRunGAMS;
    
    if iRunGAMS==1
        %Calculate the number of solutions the selected fixed and non-fixed
        %variables and ranges can generate
        dNumSols = ((1-vFixedActions).*((vMaxResult - vMinResult)./vStep+1));
        %convert zeros to ones
        for i=1:n-nObjs
            if dNumSols(i)==0
                dNumSols(i)=1;
            end
        end
        dNumSols = prod(dNumSols);

        %write to the file
        fInt = fopen('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\matlab.inc','w');
        fprintf(fInt,'*  Values set by user in MatLab\n');
        fprintf(fInt,'   SET sol Near Optimal Solutions /sol1*sol%.0f/;\n', dNumSols);
        fprintf(fInt,'   GAMMA = %s;\n',ToleranceValue);
        fprintf(fInt,'   FOPT = %.2f;\n', mMatrix(1,iOpt)/1.41);
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

        fclose(fInt);
    
    
        %now run GAMS
        %!"C:\program files\gams22.7\gams.exe"
        %!"C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\UtilityTwoStageMYNeOp.gms"
    
        [status,result] = dos('C:\Rosenberg\PhD_projects\JordanRiverBasin\UtilityIWRM\Modeling\rumyneop.bat &','-echo');
        h = msgbox('Press OK once GAMS finishes running','Running GAMS');
        uiwait(h);
        fprintf('Started next');
    end
    
    LoadJordanUtilityDataTxt;
    
    %add objective function row to vFixedActions and vFixedValues
    
    size(vFixed);
    size(vFixedActions);
    nObjs;
    
    vFixedActions = [vFixed(1:nObjs)'; vFixedActions];
    [mF nF] = size(vFixedVals);
    [mFV nFV] = size(vFixedValues);
    if mF > nF
        vFixedValues = [vFixedVals(1:nObjs); vFixedValues];
    else
        vFixedValues = [vFixedVals(1:nObjs)'; vFixedValues];
    end
    %vActLong = ['f(X)' ; vActLong];
    
    size(mObjs);
    size(mActVolumes);
    size(vFixedActions);
    size(vFixedValues);
    size(vXSteps');
    
    varargout = aSetField(varargin,'tolerance',Tolerance(2),'vFixed',vFixedActions,'vFixedVals',vFixedValues,'vStep',vXSteps','vObjLabels',vObjLabels,'vXLabels',vActLong,'vXLabelsShort',vActAbrev,'yAxisLabels',yAxisLabels);
    nearoptplotmo2(mObjs(:,1:nObjs), mActVolumes,varargout{:});
    %nearoptplotmo(hWind, mObjs(:,1:nObjs), mActVolumes, Tolerance(2), vFixedActions, vFixedValues, vXSteps', vObjLabels, vActLong, vActAbrev, yAxisLabels, fontsize, iOpt, mActCat, vGroup, 1);
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

function LoadPareto(hind,event,mMatrix,nObjs,txtTolerance,cbChecks,sSliders,mConvert,ckAdd,txtSubSpaceError,varargin)
%function LoadPareto(hind,event,hWind,mMatrix,nObjs,txtTolerance,cbChecks,vFixedVals,vStep,sSliders,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt, mActCat,ckAdd)
    %Loads and replots pareto optimal solutions
    %If ckAdd is checked, adds pareto solutions on top of the existing
    %solutions in mMatrix
    %If ckAdd is unchecked, only plots pareto solutions
    %graph
    
    [mIn nIn] = size(mMatrix);
        
    %Collect the values from the controls
    ToleranceValue = str2num(get(txtTolerance,'String'));
    SubspaceErrorNew = GetCheckTxtValue('LoadPareto',txtSubSpaceError);
    
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
    varargout = aSetField(varargin,'tolerance',ToleranceValue,'vFixed',vFixed,'vFixedVals',vFixedVals,'vStep',vStep, 'GroupToHighlight',iOpt,'SubSpaceGroup',SubSpaceGroup,'vGroup',vGroupNew, 'SubSapceError', SubSpaceErrorNew);
    nearoptplotmo2(hWind, mObjs, mNewMatrix, varargout{:});
end

function callHideSliders(hind,event,sSlider)
    %hides all sliders when checked, makes them visible when unchecked
    
    sSlider;
    
    sValue = get(hind,'Value');
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

function HideCheckboxes(hind,event,cbChecks,hTexts,hTabContainer,yPosUnchecked,yPosChecked,yPosGroupOffset, yPosGroupOffsetUnits)
    %when checked, hides all the checkboxes below the axes and moves the
    %axis text labels up to the yPosChecked value.
    
    %When unchecked, makes the checkboxes visible and returns the text
    %labels to the yPosUnchecked value
    
    % yPosGroupOffset -- the measure in yPostGroupOffsetUnits by which to
    % vertically shift the table of group text labels
      
    sValue = get(hind,'Value');
    
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
    %nG = max(size(hTextGroup));
    %nL= max(size(hLines));
    
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

%function TestAllRegion(hind,event,hWind,mMatrix,nObjs,txtTolerance,cbChecks,vFixedVals,vStep,sSliders,mConvert,vObjLabels,vXLabels,vXLabelsShort,yAxisLabels,fontsize,iOpt,mActCat,vGroup,iRunGAMS)
%Resample(hWind,event,hWindCurr,hControls,mConvert,mData,nO)
function TestAllRegion(hind,event,hWind,mMatrix,nObjs,hControls,mConvert,iRunGAMS)
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
    [txtTolerance txtNumSamples, cbChecks, txtSubSpaceError, sSliders, txtGroupOrders,txtGroupNames,cbGroupChecks] = aGetField(hControls,{'Tolerance','NumSamples','AxisChecked','SubSpaceError','Sliders','GroupOrders','GroupNames','GroupChecks'});

    ToleranceValue = get(txtTolerance,'String');  
    NumSamples = GetCheckSubSpaceError(txtNumSamples,'TestAllRegion');
    
    %size(cbChecks)
    vFixed = ReadCheckboxValues(cbChecks);
    [m n] = size(vFixed);
    
    %Calc number of decision variables
    nD = n-nObjs;
        
    %Read in the variable arguments
    varargin_fresh = getappdata(hWind,'varargs');
    [vFixedVals vGroup SubSpaceGroup TickSize iOpt vStep vXLabelsShort] = aGetField(varargin_fresh,{'vFixedVals' 'vGroup' 'SubSpaceGroup' 'TickSize' 'GroupToHighlight' 'vStep' 'vXLabelsShort'});   

    %Collect the values from the controls
    vFixedActions = vFixed((nObjs+1):n);
    
    [mObjs, mData] = SplitMatrix(mMatrix,nObjs);

    vMaxResult = max(mData);
    vMinResult = min(mData); %'C:\Rosenberg\USU\Research\CAREER\NearOptimal\EchoBMP-LP
    
    [vFixedVals vSliderSteps vFixedValsPU] = ReadSliderValues(sSliders,vFixedVals,mConvert);
    
    vStep;
    
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
        [dSols iSols] = sort(dNumSols)
        
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
     
    %nearoptplotmo(hWind, mFinalObjs, mFinalResult, ToleranceValue, vFixed, vFixedVals, vStep, vObjLabels, vXLabels, vXLabelsShort, yAxisLabels, fontsize, iOpt, mActCat, vGroupNew);
    %varargin_fresh{:}
    %vGroupNew
    varargout = aSetField(varargin_fresh,'tolerance',ToleranceValue,'vFixed',vFixed,'vFixedVals',vFixedVals,'vGroup',vGroupNew,'mGroupData',mGroupDataNewSort(:,1:2),'NumSamples',NumSamples);
    nearoptplotmo2(mFinalObjs, mFinalResult, varargout{:});
end



    