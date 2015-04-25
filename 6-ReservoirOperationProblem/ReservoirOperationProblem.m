function [hReturn vParams] = ReservoirOperationProblem(ProblemData)
% ReservoirOperationProblem.m
%
% Finds the optimal solution to a hydropower generation and irrigation reservoir operation problem
% and plots the results in the parallel coordinate plotting tool for near-optimal analysis. The optimal solution
% maximizes total (hydropower and irrigation) benefits subject to constraints on reservoir mass balance,
% reservoir capacity, hydropower turbine capacity, river mass balance, and maintaining an instream flow at a location downstream of
% the irrigation diversion. Hydropower turbine releases may also be used for irrigation.
% User provides input data describing reservoir inflows, initial reservoir storage, reservoir and turbine capacities, and
% minimum required flows at river location A. A schematic for the problem:
%
%                         Spill
%  Inflow           |----------\
%  ---------->      |-Turbine--------------> A --->
%        Reservoir  |                 \
%                                      \ Irrigation

% INPUT PARAMETERS
%  ProblemData = a structure with the problem data in the fields (empty or
%                missing mean the default value is used)
    %  inflows = 1 x T row vector of reservoir inflow in arbitrary volume
    %       units (e.g., Mm3) (T=the number of time steps) (Deafult: [2 2 3 4 3
    %       2])
    %  initStorage = Initial reservoir storage in arbirtrary volume units
    %       (e.g., Mm3) (Default: 5)
    %  maxStorage = Reservoir capacity in arbirtrary volume units (e.g., Mm3)
    %       (Default: 10)
    %  minFlowAtA = Minimum required flow at river location A in arbitratry
    %       flow units (e.g., Mm3/time step) (Default: 1)
    %  TurbineCapacity = Maximum release through the hydropower turbine in
    %       arbitrary flow units (e.g., Mm3/time step) (Default: 4)
    %  hydropowerBenefits = 1 x T row vector of Hydropower benefits in
    %       arbitrary money/volume units (e.g., $/Mm3) (Default value: [1.6 1.7
    %       1.8 1.9 2.0 2.0])
    %  irrigationBenefits = 1 x T row vector of Irrigation benefits in arbitrary money/volume
    %       units (e.g., $/Mm3) (Default value: [1.0 1.2 1.9 2.0 2.2 2.2])
%    
% OUTPUTS
%   hReturn = handle to the parallel coordinate plot generated       
%   vParams = cell array of parameter value pairs passed to the plot function
%
%  CALLED FUNCTIONS
%    - toolbox :  \optim\optim\linprog.m
%    - stratgibbs.m (stratify Gibbs sample from polytope)
%    - nearoptplotmo2.m (plot results in interactive parallel coordinate
%       plot tool


% To run, download all the files in the folders 2-GenerateAlternatives and
% 3-InteractiveParallelPlot on the GitHub repository.
%
% #####################
%   Programmed by David E. Rosenberg
%   April 2015
%   
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   Citation:
%   David E. Rosenberg (2015) "Blended near-optimal alternative generation, 
%   visualization, and interaction for water resources decision making".
%   Water Resources Research. doi:10.1002/2013WR014667.
%   http://onlinelibrary.wiley.com/doi/10.1002/2013WR014667/full

%   Licensing:
%   This code is distributed AS-IS with no expressed or implied warranty regarding the claimed functionality. The entire code or parts 
%   may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires 
%   prior written permission from the author.
%
%   Bug Reports and Feedback:
%   Bug reports and Feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note that while much appreciated, there is no promise of when, or if, a reported bug will be corrected.


%% Step A. Define the Linear Program.
%
% Our problem is Max Z = cX
%                st Aeq*X = Beq
%                   X <= uB
%                   X >= lB

   %Default return values
   mDecs = []; mObjs=[]; vParams=[];

   %Define the defalt problem data
   ToleranceValue = 1; %Dummy placeholder for near-optimal tolerance value
   dfProblemData.inflows = [2 2 3 4 3 2]; % reservoir inflow in arbitrary volume units (e.g., Mm3)
   dfProblemData.hydropowerBenefits = [1.6 1.7 1.8 1.9 2.0 2.0]; % Hydropower benefits in arbitrary money/volume units (e.g., $/Mm3)
   dfProblemData.irrigationBenefits = [1.0 1.2 1.9 2.0 2.2 2.2]; % Irrigation benefits in arbitrary money/volume units (e.g., $/Mm3)
   dfProblemData.initStorage = 5; % Initial reservoir storage in arbirtrary volume units (e.g., Mm3)
   dfProblemData.maxStorage = 10; % Reservoir capacity in arbirtrary volume units (e.g., Mm3)
   dfProblemData.minFlowAtA = 1; % Minimum required flow at river location A in arbitratry volume units (e.g., Mm3)
   dfProblemData.TurbineCapacity = 4; % Maximum release through the hydropower turbine 
   
   %Build the actual problem data from data provided (and default values,
   %    when field(s) are missing
   if nargin < 1 || ~isstruct(ProblemData)
       ProblemData = dfProblemData;
   else
       %Check each field of the structure. If any missing, read in the
       %default value
       if ~isfield(ProblemData,'inflows')
           ProblemData.inflows = dfProblemData.inflows;
       end
       
       if ~isfield(ProblemData,'hydropowerBenefits')
           ProblemData.hydropowerBenefits = dfProblemData.hydropowerBenefits;
       end
       
       if ~isfield(ProblemData,'irrigationBenefits')
           ProblemData.irrigationBenefits = dfProblemData.irrigationBenefits;
       end
       
       if ~isfield(ProblemData,'initStorage')
           ProblemData.initStorage = dfProblemData.initStorage;
       end
  
       if ~isfield(ProblemData,'maxStorage')
           ProblemData.maxStorage = dfProblemData.maxStorage;
       end
       
       if ~isfield(ProblemData,'minFlowAtA')
           ProblemData.minFlowAtA = dfProblemData.minFlowAtA;
       end
       
       if ~isfield(ProblemData,'TurbineCapacity')
           ProblemData.TurbineCapacity = dfProblemData.TurbineCapacity;
       end
   end

   %final check on integretedy of the data structure
   T = size(ProblemData.inflows,2);
   if (T~=size(ProblemData.irrigationBenefits,2)) || (T~=size(ProblemData.hydropowerBenefits,2))
       error('Input data of inconsistent size. Check inflows, irrigationBenefits, and/or hydropowerBenefits')
       
       return
   end
   
   %Decision Variables
   Stor = zeros(1,T); %Reservoir storage
   Spill = zeros(1,T); %Reservoir spill
   HydroRel = zeros(1,T); %Reservoir release to generate hydropower
   Irrigat = zeros(1,T); % Delivery to irrigation area
   FlowAtA = zeros(1,T); % Instream flow at point A on the river
   
   X = [Stor Spill HydroRel Irrigat FlowAtA]; %Concatenation of all variables
   c = [zeros(1,2*T) ProblemData.hydropowerBenefits ProblemData.irrigationBenefits zeros(1,T)]; %Concatenation of objective function coefficients
                % that sum-totals hydropower and irrigation benefits
                
   
   
   %Constraints and bounds
   %  - Equity constraints
   %     - Reservoir storage balance in time step 1
         Aresbalinit = [1 zeros(1,T-1) 1 zeros(1,T-1) 1 zeros(1,T-1+2*T)];
   %     - Reservoir storage balance in time steps 2 to T
         Aresbal = [-1 1 zeros(1,T-1) 1 zeros(1,T-1) 1 zeros(1,T-2+2*T)];
   %     - River balance in each time step (rows T+1 to 2T)
         Arivbal = [zeros(1,T) 1 zeros(1,T-1) 1 zeros(1,T-1) -1 zeros(1,T-1) -1 zeros(1,T-1)];
         
   %Combine into a single matrix so the constraints are of the 
   %form Aeq*X = Beq
   
   Aeq = Aresbalinit;
   for i=1:T-1
       Aeq = [Aeq; circshift(Aresbal',i-1)'];
   end
   for i=1:T
       Aeq = [Aeq; circshift(Arivbal',i-1)'];
   end
           
   Beq = [ProblemData.inflows zeros(1,T)]';
   Beq(1) = Beq(1)+ProblemData.initStorage;
   
   %Bounds on decision variables
   %Lower bounds
   lB = [zeros(T*4,1); (ProblemData.minFlowAtA)*ones(T,1)];
   %ending reservoir storage must equal or exceed initial storage
   lB(T) = ProblemData.initStorage;
   
   %Upper bounds
   uB = [(ProblemData.maxStorage)*ones(1,T) inf(1,T) (ProblemData.TurbineCapacity)*ones(1,T) inf(1,2*T)]';
   
   n=size(X,2); %Number of decision variables
   m=size(Aeq,1); %Number of constraints
   

%% Step B. Solve for the optimal solution
    [Xopt,fopt] = linprog(-c,[],[],Aeq,Beq,lB,uB); % negative sign on c indicates maximization   

    sprintf('Optimal objective function value: %.2f\nOptimal solution:',-fopt)
    Xopt

%% Step C. Plot optimal solution in parallel coordinate tool and allow for interactive
% near-optimal exploration

%Generate the Near-optimal tolerance constraint (cX >= f* ToleranceValue)
%for a maximization problem
    A_nearopt = [-c];%[-eye(n);-c]; %Add dummy non-negativity constraints.
    b_nearopt = [fopt*ToleranceValue];%[zeros(n,1);fopt*ToleranceValue];

    %Put the model into ProblemStruct
    ProbStruct.Aineq = A_nearopt;
    ProbStruct.bineq = b_nearopt;
    ProbStruct.Aeq = Aeq;
    ProbStruct.beq = Beq;
    ProbStruct.ub = uB;
    ProbStruct.lb = lB;


    % Plot the sampled alternatives and optimal solution using the interactive parallel-cordinate plotting tool nearoptplotmo2.m and allow the user to
    % interact. Notes:
    %    - Optimal solution in row 1 of the data sets for the objective function and decision varaibles, Sampled alternatives in rows 2 to p+1
    %    - Reverse sign on optimal objective function value because it's a maximization problem
    %    - Add the near-optimal tolerance constraint as row 1 of the inequality constraint set. This constraint is of the form:
    %        c*X >= Zopt*Tolerance
    %      which is transformed into
    %        -c*X <= -Zopt*Tolerance
    %    - Render traces in same colors across objective and decision variable
    %         axes
    %    - Organize the groups so they plot sampled alternatives on bottom and optimal solution last
    %       on top and highlighted in thick black.
    %

    vGroups = {'Optimum'};
    mGroupData = ['Optimum' {1} {2.5}];
    mActCat = [repmat({'Storage'},T,1); repmat({'Spill'},T,1); repmat({'Hydro Release'},T,1); ...
           repmat({'Irrigation'},T,1) ; repmat({'Flow at A'},T,1)];
    %Build the time labels
    vTimes = cell(1,T);
    for i=1:T
        vTimes{i} = sprintf('%s %d','Time Step',i);
    end
    vXLabels = repmat(vTimes,1,5);
    vObjLabels = {'Total Benefits'};    
    yAxisLabels = {'Benefits (arbitrary units)', 'Water volume (arbitrary units)'};

    %Build the maxtrices of solutions/alternatives
    vParams = {'Tolerance',ToleranceValue,'ProbForm',ProbStruct,'cFunc',c,...
          'OptSolRow',1,'NearOptConstraint',size(A_nearopt,1),'mActCat',mActCat,'ShowControls',0,'FontSize',20,'YAxisMargins',[0.25 0.25],...
          'vXLabels',vXLabels,'vObjLabels',vObjLabels,'yAxisLabels',yAxisLabels,'HideCheckboxes',1, ...
          'ShowObjsDiffColor',0,'vGroup',vGroups,'mGroupData',mGroupData,'GroupToHighlight','Optimum', ...
          'ShowGroupLabels',1,'ErrorResid',-1e-6,'StartTab',1,'GenerateType',3,'GenerateMethod',2, ...
          'AxisScales','custom',[1 1],[0 zeros(1,T*5);60*ones(1,T*5+1)]};

    %The plot command
    nearoptplotmo2(-fopt,Xopt',vParams{:});

end

