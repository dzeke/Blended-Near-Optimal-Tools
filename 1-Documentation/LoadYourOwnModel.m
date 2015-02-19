% LoadYourOwnModel.m
%
% Provides an example of how to load your own linear optimization program into the
% near-optimal tools and use the tools. Also
% provides guidance for other types of models like mixed-integer and
% non-linear programs.
%
% To run, download all the files in the folders 2-GenerateAlternatives and
% 3-InteractiveParallelPlot on the GitHub repository.
%
% #####################
%   Programmed by David E. Rosenberg
%   August 2014
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
%   may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires 
%   prior written permission from the author.
%
%   Bug Reports and Feedback:
%   Bug reports and Feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note that while much appreciated, there is no promise of when, or if, a reported bug will be corrected.


%% Loading an example Linear Program.
% 
%  This example linear program is an n=5-decision variable problem (x1, x2, ..., x5) to maximize the linear-combination of the
%  variables (Z = 5*x1 + 3*x2 + x3 + 0.5*x4 + 0.5*x5) subject to constraint
%  that solutions fall within a 5-dimensional simplex with side length r
%  (i.e., X1, X2, ..., Xn >= 0, X1 + X2 + ... + X5 <= r) 
%
%  The model is 
%    maximize     Z = cX
%    such that    Ax <= b
%  Where
%   c = 1 x n row vector of objective function coeffiecients
%   A = m x n matrix of contraint coeffients that define the simplex
%   b = m x 1 column vector of right-hand side coefficients of constraints that define the simplex
%   n = number of decision variables (length of c and b and rank of A)
%   r = length of a simplex side
%   p = number of alternatives to sample
%
%
% The Basic Steps are:
%  A. Define the model and model data
%  B. Solve for the optimal solution and sample alternatives
%  C. Plot them up and allow the user to interact

% Step A. Defining the model and model data
n = 5;                  %Number of decision variables
c = [5 3 1 0.5 0.5];    %
r = 6;
A = [-eye(n); ones(1,n)]; %first n rows represent the greater than or equal to inequalities ; (negative sign reverses the direction of default less than)
                          %last row is X1 + X2 + ... + Xn <= r
b = [zeros(n,1); r];    % right hand side of constraints
p = 150;       % Number of alternatives to generate

mVert = [zeros(1,n); r*eye(n)]; %vertices (corner points) of the simplex. Used for comparison on the plot.

% Step B. Solve for the optimal solution and sample alternatives
[Xopt,fopt] = linprog(-c,A,b); % negative sign on c indicates maximization

%Generate the Near-optimal tolerance constraint (cX >= f* ToleranceValue)
%for an maximization problem
ToleranceVal = 0;
A_nearopt = [A;-c];
b_nearopt = [b;fopt*ToleranceVal];

%Put the model into ProblemStruct
ProbStruct.Aineq = A_nearopt;
ProbStruct.bineq = b_nearopt;

% Random sample alternatives from within the feasible region using the
% stratified gibbs tool. Notes:
%    - This sampling approach works for any closed,bounded region, not just near-optimal!
%    - Will also stratify along specified linear combinations of decision variables --  in this case the objective function c
%    - Use the optimal solution Xopt as a starting point
[Xs,vVal] = stratgibbs(p,ProbStruct,struct('lincombo',c','x0',Xopt','extmethod','opt','errorresid',0));

% Only use valid alternatives (vVal>0)
XsValid = Xs(vVal>0,:);

% Calculate objective function values for the valid sampled alternatives. 
vObjs = XsValid*c';
% Calculate objective function values for the vertices
vObjsVert = mVert*c';

% Step C. Plot the sampled alternatives, vertices, and optimal solution up using the interactive parallel-cordinate plotting tool nearoptplotmo2.m and allow the user to
% interact. Notes:
%    - Optimal solution in row 1 of the data sets for the objective function and decision varaibles, Sampled alternatives in rows 2 to p+1, 
%       and Vertices in the final n+1 rows.
%    - Reverse sign on optimal objective function value because it's a maximization problem
%    - Consider a near-optimal tolerance of 85% of the optimal objective function value
%    - Add the near-optimal tolerance constraint as row 7 of the constraint set. This constraint is of the form:
%        c*X >= Zopt*Tolerance
%      which is transformed into
%        -c*X <= -Zopt*Tolerance
%    - Render traces in same colors across objective and decision variable
%         axes
%    - Organize the groups so they plot sampled alternatives on bottom, vertices next with a thicker line, and optimal solution last
%       on top and highlighted in thick black.
%

vGroups = ['Optimum'; repmat({'Feasible region'},length(vObjs),1); repmat({'Vertices'},length(vObjsVert),1)];
mGroupData = ['Feasible region' {1} {1}; 'Vertices' {1} {1.5}; 'Optimum' {1} {2.5}];
nearoptplotmo2([-fopt;vObjs;vObjsVert],[Xopt';XsValid;mVert],'Tolerance',ToleranceVal,'ProbForm',ProbStruct,'cFunc',c,...
      'OptSolRow',1,'NearOptConstraint',size(A_nearopt,1),'ShowControls',0,'FontSize',20,'YAxisMargins',[0.25 0.25],...
      'ShowObjsDiffColor',0,'vGroup',vGroups,'mGroupData',mGroupData,'GroupToHighlight','Optimum');

% Suggested interactions to emphasize some key aspects of the example:
% 1. To show near-optimal alternatives, select the interact tab. In the
%       Generate New Alternatives box, set type to Random sample. The Engine will
%       update to MATLAB (LP matrix). Enter the # Samples (e.g., 250) and click
%       Generate. New near-optimal alterantives will appear in purple on the
%       plot and show the portion of the feasible region (simplex) where
%       near-optimal alternatives reside.
%     - Click the Display tab and rename the second group (Feasible
%     region_Resample) to Near-optimal.
% 2. To show how the near-optimal alternatives vary with the objective
% function value,
%     - On the Display Tab, uncheck the first and third 
%       groups (Feasible region and Vertices) and click Show Checked Groups. Only
%       traces for the Near-optimal and opitmum groups will show.
%     - On the plot, check the box above the label Objective 1
%     - Click the Color Ramp tab, enter the # Color classes as 5, keep
%       Direction at Ascend, and check the box Ramp color. The plot will update and
%       show the darkest green lines have large values for Decision variable
%       1 (direct relationship), varied but generally small values for Decision
%       variable 2, and near-zero values for the remaining decision
%       variables.

%% Loading results for Mixed-Integer problems

% For mixed-integer problems, load results for the optimal solution and
% sampled alternatives as in step C above for the linear program with the following changes:
%   a. Exclude the parameters AMat, Brhs, cFunc, OptSolRow, and NearOptConstraint. 
%   b. Specify the GAMS file where the mixed-integer model resides
%   (myGamsFile)
%   c. Specify the list of set element labels for the decision variables in the
%   GAMS file. That is, these labels map each column in mActs (a decision
%   variable) to a decision variable in the GAMS model.
%   d. Specify the Generate Type as 4 (enumerate) and Generate Method as 3 (GAMS) to set the
%   Generate Type and Engine fields correctly in the Generate New Alternatives box
%   on the Interact Tab when the plot loads.

% 
%     >> nearoptplotmo2(mObjs,mActs,'Tolerance',Gamma,'GroupToHighlight',vGroupText{1}, ...
%                'vGroup',vGroupText,'sGamsFile',myGamsFile, 'vXLabelsShort',vActNew, ...
%                'GenerateType',4,'GenerateMethod',3);

% To correctly pass data to GAMS, your GAMS file
% must read the following GAMS inputs from the file matdata2.gdx:
%   - GAMMA (near-optimal tolernace)
%   - LFixed (nD x 1)
%   - LFixedVal
%   - LFixedValVol
%   - InputPassed
%   - RunMode
%
% To correctly pass data from GAMS back to Matlab, your GAMS file must
% write the following GAMS outputs to the file matsol.gdx:
%  - Costs
%  - LongActs
%  - LevelData
%  - LToUse
%  - NumSols
%  - NumSolvs
%  - vTurnedFixed

% See further explanations of the use of these GAMS parameters in the
% Matlab wrapper \3-InteractiveParallelPlot\EnumNEIntSolsGams4.m. See also 
% a more detailed example in folder E - in the application for Amman,
% Jordan.

%% Loading results for non-linear programs

% Currently, the tools only support the interactive plotting and rendering
% of pre-calculated results for non-linear programs. This is an active area
% of further research. Use the current tools as follows:

% >> nearoptplotmo2(mObjNonLinear,mDecsNonLinear)

% where mObjNonLinear and mDecsNonLinear are vectors/matrices whose rows
% define non-linear solutions and alternatives. Only tools on the Plot Data
% and Controls menus as well as the Color Ramp and Display tabs will work.