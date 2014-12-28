function [xMGA, p, MinDist, mTrialInfo, StopCode, StopExplain, runTime] = doMGA(xStart,options,optopoptions)
% Generates near-optimal alternatives for a base optimization problem using the specified Modeling to Generate
% Alternatives (MGA) method.
%
%  The supported MGA methods are:
%     1 = Hop-Skip-Jump that searches for a maximally different
%               alternative by minimizing the sum of the non-zero variables
%               in all prior alternatives and continues until no new variables enter
%               the basis or all the decision variables are in the basis.
%               See Brill, E. D., Jr., Chang, S.-Y., and Hopkins, L. D. (1982).
%               "Modeling to Generate Alternatives: The HSJ Approach and an Illustration Using a Problem in Land Use Planning." Management Science, 28(3), 221-235.
%
%     2 = MGA-serial -- Distance norm that searches for a maximally different
%               alternative by maximizing the absolute value of the distance to the nearest prior
%               alternative (L-infinity norm). Run in an iterative approach until a stopping criterial is reached (run time, minimum distance, number of alternatives).
%               See Zechman, E. M., and R. Ranjithan (2007), Generating Alternatives Using Evolutionary Algorithms 
%               for Water Resources and Environmental Management Problems, Journal of Water Resources Planning and Management, 133(2), 156-165. 
%               http://ascelibrary.org/doi/abs/10.1061/%28ASCE%290733-9496%282007%29133%3A2%28156%29.
%
%     3 = MGA-simultaneous -- Distance norm that searches for maximally different alternatives by maximizing the absolute value of the
%               distance to the nearest neighboring alternative (L-infinity norm). Same as #2 except run simultaneously (rather than iteratively) 
%               so that all alternatives identified at once (much larger program than #2).
%               See Zechman, E. M., and R. Ranjithan (2007), Generating Alternatives Using Evolutionary Algorithms 
%               for Water Resources and Environmental Management Problems, Journal of Water Resources Planning and Management, 133(2), 156-165. 
%               http://ascelibrary.org/doi/abs/10.1061/%28ASCE%290733-9496%282007%29133%3A2%28156%29.
%
%     Note slight modification here for #2 and #3 as we include optimal
%     solution in the set over which distances are calculated since we have
%     already identified the optimal solution.
%
%   The base optimization problem is
%       Max (or Min) objfun(x)
%       s.t. A x <= b
%            Aeq x = beq
%            nlconfun(x) <= 0 and nlconfun(x) = 0
%            x <= uB
%            x >= lB
%
%       If the objfun and tolerance are specified, doMGA will add a
%       constraint of the form 
%
%       s.t. objfun(x) <= tolerance*objfun(xOpt)
%           for a minimization problem (tolerance > 1), or
%       s.t. objfun(x) >= tolerance*objfun(xOpt)
%           for a maximization problem (0 < tolerance < 1)
%
% INPUTS
%   - xStart = a 1 x n row vector of decision variable values that define the start point (typically the optimal
%       solution to the underlying optimization problem) for alternative generation.
%
%   - options = a structure with the following fields that define both the
%       optimization problem and how to generate near-optimal alternatives with
%       MGA. The option fields are:
%
%       objfunc = objective function used with the tolerance parameter (see below) to calculate and add the near-optimal tolerance constraint to the model formulation. If omitted
%           no constraint is added
%           
%           for a linear program - an n x 1 row vector of decision variable coefficients
%           for a nonlinear program - a function handle (@...) specifying the
%               the nonlinear objective function as a function of f(x). This function must comply with
%               requriements for the objective function listed in fmincon.m. 
%
%       tolerance = real number specifying the near optimal tolerance. Used
%           with the objfun to calculate and add a near-optimal tolerable
%           deviation constraint to the model formulation
%              - Maximization problem: 0 <= tolerance <= 1 
%              - Minimization problem: tolerance > 1

%       A = m x n matrix of coefficients representing the m linear
%           inequality constraints to the problem in the form A x <= b (default value
%           is [])
%
%       b = an m x 1 column vector of the right-hand side values of the m
%           linear inequality constraints, again A x <= b (default value is [])
%
%       Aeq = q x n matrix of coefficients representing the q linear
%           inequality constraints to the problem in the form Aeq x = beq (default value
%           is [])
%
%       beq = an q x 1 column vector of the right-hand side values of the
%           linear equity constraints, again Aeq x = beq (default value is [])
%
%       nonlincon =  a function handle to the non-linear inequality and
%           non-linear equaliity constraint functions (e.g., @...). This function must comply with
%           requriements for the constraint function listed in fmincon.m.
%           (Default value: '')
%
%       lB = n x 1 vector of lower bounds for the n decision variables
%           (default value is [])
%
%       uB = n x 1 vector of upper bounds for the n decision variables
%       (default value is [])
%
%       errorcrit = scalar value determining the error criteria to use to
%                   classify whether a decision variable is non-zero in the
%                   MGAType = 1 method (iterative).
%                   i.e., abs(x(i)) > errorcrit is non-zero. Default value
%                   (0)
%
%       MGAType = a switch variable that takes the values 1, 2, or , 3 to indicate the method (default MGAType=1):
%           1 = Hop Skip Jump
%           2 = Serial
%           3 = Simultaneous  
%
%       MaxAlts = an integer to specify the maximum number of alternatives
%           to geenerate (default value: 1)
%
%       StopDistance = a single value expressed in units of the decision
%           variables that represents the distance below which to stop
%           generating alternatives
%
%       StopTime = a single value expressed in seconds that specifies the
%           maximum length of time to run MGA interations (default value: 20).
%
%       SolveAsGA = Binary that takes the value of 1 to solve the model with a Genetic Algorithm
%           (Matlab's ga). Default value: 0 (regular LP or NLP)   
%
%       IgnoreStartPoint = binary that takes the value of 1 to ignore the
%           xStart point provided (Only valid for Simultaneous method)
%           (Default: 0)               0
%
%    optoptions = a structure of optimization program option settings
%    specified using optiset for the problem.
%
%
%  OUTPUTS
%    xMGA = a p x n matrix representing the p near-optimal alternatives
%       generated.
%    p = the number of alternatives generated
%    MinDist = minimum distance between alternatives
%    mTrialInfo = p x 2 matrix of minimum distance of the generated alternative to others and cummulative runtime (for HSJ and Serial MGATypes) 
%    StopCode = an integer code representing the reason interactions were stopped
%       -2 = optimization error on solve attempt
%       -1 = error with input
%        0 = other unspecified reason
%        1 = Maximum number of alternatives reached (p==NumAlts)
%        2 = Stop Distance reached (current distance < StopDistance)
%        3 = Stop time reached (run time > StopTime)
%        4 = # basic variables equals total number of variables (HSJ
%               method)
%        5 = # basic variables on current iteration equals number of basic
%              variables on prior iteration (HSJ method)
%     StopExplain = a text description of the StopCode
%     runTime = the run time of iterations in seconds measured using tic
%       toc
%
%  CALLED FUNCTIONS
%    - toolbox :  \optim\optim\fmincon.m
%    -    ga.m
%    -    gs.m
%    -    maxextentind.m
%    -    maxextentindnl.m
% 
%% #####################
%   Programmed by David E. Rosenberg
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%   - November, 2014. First version. To compare MGA run times to new
%   near-optimal alternative generation methods.
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
%   Bug Reports and Feedback:
%   This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.
%
%% #######################

    StopCode = -1;
    
    StopCodesText = {'Optimization error' ...  % -2
        'Error with input'  ...             % -1
        'Other unspecified reason'  ...     % 0
        'At max. # of alts', ... % 1
        'At min distance'  ...    % 2
        'At max time'  ...        % 3
        '# basic variables = n'  ... %4
        'No new basic variables added'}; %5

    %Read / error check input
    if isempty(xStart)
        error('Need to define xStart (even if ignoring)')
        return
    end
    
    [r,n] = size(xStart);
    
    %Set default options
    A = []; b=[]; Aeq = []; beq = []; MGAType = 1; MaxAlts = 1; StopDistance = 0.1*max(abs((xStart(1,:))));
    StopTime = 20; lB=[]; uB=[]; nonlincon = ''; errorcrit = 0; SolveAsGA = 0;
    IgnoreStartPoint = 0;
    
    %Flags for storing state of input parameters and subsequent actions needed to take
    blAddNearOptimal = 0;
    blAsLP = 1; %As linear program 1=yes, 0=no
    blLinObjFunc = 1; %Linear objective function 1=yes, 0=no
    
    %Read in options
    if isstruct(options)
        %A and b
        if isfield(options,'A') && isfield(options,'b')
            Atmp = options.A;
            Btmp = options.b;
            
            if (size(Atmp,1) == size(Btmp,1)) && (size(Atmp,2) == n)
                A = Atmp;
                b = Btmp;
            else
                warning('options.A and options.b improperly/inconsistently sized. Ignoring. Using defaults of [].')
            end
        end
        
        %Aeq and beq
        if isfield(options,'Aeq') && isfield(options,'beq')
            Aeqtmp = options.Aeq;
            Beqtmp = options.beq;
            
            if (size(Aeqtmp,1) == size(Beqtmp,1)) && (size(Aeqtmp,2) == n)
                Aeq = Aeqtmp;
                beq = Beqtmp;
            else
                warning('options.Aeq and options.beq improperly/inconsistently sized. Ignoring. Using defaults of [].')
            end
        end   
        
        %Lower and upper bounds        
        if isfield(options,'lB')
            lBtmp = options.lB;
                     
            if (size(lBtmp,2) == n)
                lB = lBtmp;
            else
                warning('options.lB improperly/inconsistently sized. Ignoring. Using default of [].')
            end
        end   
        
        if isfield(options,'uB')
            uBtmp = options.uB;
                     
            if (size(uBtmp,2) == n)
                uB = uBtmp;
            else
                warning('options.uB improperly/inconsistently sized. Ignoring. Using default of [].')
            end
        end   
        
        if isfield(options,'nonlincon')
            nonlincon = options.nonlincon;
            blAsLP = 0;
        end
        
        %MGA iteration criteria
        if isfield(options,'MGAType') 
            if ismember(options.MGAType, [1 2 3])
                MGAType = options.MGAType;
            else
                warning('options.MGAType outside acceptable values. Ignoring. Using default of 1 (Hop-Skip-Jump)')
            end
        end        

        if isfield(options,'MaxAlts') 
            MaxAlts = round(options.MaxAlts);
        end          
        if isfield(options,'StopDistance') 
            StopDistance = options.StopDistance;
        end    
        if isfield(options,'StopTime') 
            StopTime = options.StopTime;
        end   
        if isfield(options,'errorcrit')
            errorcrit = options.errorcrit;
        end
        if isfield(options,'SolveAsGA')
            SolveAsGA = options.SolveAsGA;
        end
        if isfield(options,'IgnoreStartPoint')
            IgnoreStartPoint = options.IgnoreStartPoint;
            if (IgnoreStartPoint==1) && (MGAType < 3)
                warning('Must use specified start point for MGA Method %d',MGAType)
            end
        end       
        
    end
    
    if isfield(options,'objfunc') && isfield(options,'tolerance') && (options.tolerance >=0)
        % An objective function and near-optimal tolerance are specified. 
        % Expand the constraints to include the near-optimal tolerance
        % constraint
        
        objfunc = options.objfunc;
        tolerance = options.tolerance;
        
        if isnumeric(objfunc)
            %linear program
            %check size consistency
            
            if size(objfunc,1) == n     
               %add near-optimal tolerance constraint in the form for a
               %minimization problem
               A = [A; objfunc']; 
               b = [b; tolerance*xStart*objfunc];                
                                
               if tolerance < 1 
                    % Reverse direction for a maximization problem
                    A(size(A,1),:) = -A(size(A,1),:);
                    b(size(b,1),:) = -b(size(b,1),:);
               end
               blAddNearOptimal = 0;
            else
                warning('objfun is not properly sized. Continuing with default of no near-optimal constraint')
                objfunc = [];
            end
        elseif strcmpi('function_handle',class(objfunc))
            %Nonlinear formulation
            blAsLP = 0;
            blLinObjFunc = 0;
            NEval = tolerance*objfunc(xStart);
            
            %Condition cFunc and NEval based on the direction of the
            %optimization
            direction = 1; %minimization
            if tolerance < 1
                %maximization
                direction = -1;
            end
            
            if strcmpi('char',class(nonlincon))
                %empty, no nonlinear constraints. Just need near-optimal
                %tolerance
                nonlincon = @(x)NearOptNLConstraints(x,objfunc,NEval,direction);
            else
                %Pre existing non-linear constraints, include through
                %function handle
                nonlincon = @(x)NearOptNLConstraints(x,objfunc,NEval,direction,nonlincon);
            end
        end
    else
        if isfield(options,'objfunc') || isfield(options,'tolerance')
            warning('Need to specify both an objective function and a tolerance to include a near-optimal constraint. Continuing with default of no near-optimal constraint')
        end
    end   
    
    %Update nlconfun in options sturcture
    options.nonlincon = nonlincon;
    
    %Prepare the outputs
    xMGA = xStart(1,:);
    MinDist = 0;
    mTrialInfo = [];
    p = 0;
    StopCode = 0;   
    runTime = 0;
    
    %% Start MGA iterations
    
    runTime = 0;
    tic;
    
    if MGAType==1
      % Work by the Hop-Skip-Jump method that minimizes the sum of all
      % non-zero variables in all prior alternatives.
      NonZeroDecs = (abs(xMGA)>errorcrit)+0; %The non-zero deicision variables. Take largest value of any prior altnernative. Once a decision varuabke is non-zero, always non-zero. 0 at end to turn to numeric
      NumNonZerosDecs = sum(NonZeroDecs,2);
       
      %Continue generating new MGAs while
      %  1) number of non-zero decision variables is less than number of
      %     decision variables, or
      %  2) number of non-zero decision variables in current alternative is
      %     greater than the number of non-zero decision variables in the
      %     prior alterative 
      %  i.e., stopping criteria are not yet satisified
      
      while StopCode == 0
          pNext = p+1; %Temporary advance counter
          
          if SolveAsGA % Solve with a genetic algorithm
            options = gaoptimset('MutationFcn',@mutationadaptfeasible,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr,@gaplotrange}, 'Display','iter','TimeLimit',StopTime-runTime);  
            [xNew,fBest,errorflag,Test2] = ga(@(x)x*NonZeroDecs(pNext,:)',size(xMGA,2),A,b,Aeq,beq,lB,uB,nonlincon,options); 
            
            Test2
            
          elseif blAsLP  %Solve with a linear program
            [xNew, fNew, errorflag, Test2] = linprog(NonZeroDecs(pNext,:),A,b,Aeq,beq,lB,uB,[],optimoptions('linprog','maxiter',15000,'Display', 'off','algorithm','simplex'));
            xNew = xNew';
          else %Solve as a non-linear program
            [xNew,fNew,errorflag,Test2] = fmincon(@(x)x*NonZeroDecs(pNext,:)',xMGA(pNext,:),A,b,Aeq,beq,lB,uB,nonlincon,struct('maxiter',15000,'Display', 'off'));
          end
          
          runTime = toc;
           
          if errorflag > 0
            %Record the new alternative
            p=pNext; %Increment the alternative counter
            xMGA = [xMGA;xNew];            
            
            %Calculate distance of current alternative to prior ones
            if p==1
                mComp = [];
            else
                mComp = [1 2];
            end
            MinDist = DistanceToNearestNeighbor(reshape(xMGA(1:end-1,:),n*p,1)',p,n,mComp,xNew);
            %Log this distance for this iteration
            
            mTrialInfo = [mTrialInfo;MinDist runTime];
            NonZeroDecs = [NonZeroDecs;(max(abs(xMGA))+0>errorcrit)]; %
            NumNonZerosDecs = [NumNonZerosDecs; sum(NonZeroDecs(p+1,:),2)];
          else
              Test2.message
              StopCode = -2;
          end
          
          
          
          %Check stopping criteria
          if p>=MaxAlts
              StopCode = 1;
          elseif runTime > StopTime
              StopCode = 3;
          elseif (NumNonZerosDecs(p+1)>=n)
              StopCode = 4;
          elseif (NumNonZerosDecs(p+1) <= NumNonZerosDecs(p))
              StopCode = 5;
          end
      end 
      
      %[[1:p+1]' NumNonZerosDecs NonZeroDecs]
   
    elseif MGAType==2
      % Run by iteratively and serially searching for a maximally different
      % alternative by maximizing the absolute value of the distance to the nearest prior
      % alternative (L-infinity norm) until a stopping criterial is reached (run time, minimum distance, number of alternatives)
      %
      % The compact problem is
      %  Max Z = min(sum(abs(x_k - x_k,z))))
      %   s.t. Original linear inequality, equality, and non-linear
      %   constraints
      %    where k= decision variable index (1..n)
      %          z= index of optimal solution and prior generated
      %          alternatives (0, 1, .. p-1)
      %
      % A larger, Full version of the problem incluedes variables for the minimum distances between
      %             the current alternative and each prior alternative:
      % Max Z
      % s.t Original Linear inequality, Equality, and non-linear constraints
      %    Z <= Dist(z)   for all z <= p-1    [L-infinity norm]
      %    Dist(z) <= sum_i(abs(x(i) - VERTS(i,z)))   for all z, [distance
      %                                                             measure]
      %    Dist(z) <= sum_i(MaxExtent(i)-MinExtent(i))    for all z,  [additional upper
      %                                                 bound that is the
      %                                                 most a variable can
      %                                                 change
      %
      % z expands by one member at each iteration. 
      
      % First define the upper bound for distance variables 
      % Turn off the objfunc option before calling maxextentind/nl
      if isstruct(options) && isfield(options,'objfunc')
          optionsForExtent = rmfield(options,'objfunc');
      else
          optionsForExtent = options;
      end
      
      if blAsLP  
          [xExtsFull,xExtsCompact] = maxextentind(A,b,optionsForExtent);
          UpperBoundDist = sum(xExtsCompact(:,2)-xExtsCompact(:,1));
      else
          [xExtsFull,xExtsCompact] = maxextentindnl(A,b,optionsForExtent);
          UpperBoundDist = sum(xExtsCompact(:,2)-xExtsCompact(:,1)) ;           
      end
      
      %Formulate the augmented problem
      %Organize the decisions as:
      %   Z  x(1..i)   Dist(1..v)
      %Organize the constraints as:
      %  Ax  <= b
      %  Z <=  Dist(1..v)  
      %  Dist(v) <= UpperBoundDist
      
      [m,n] = size(A);
      if SolveAsGA
          xMGAFull = xMGA;
      else
          xMGAFull = [0 xMGA]; %first column is the minimum distance variable
      end
      
      while StopCode == 0
          %Build the matrix inputs describing the problem with pNext alternatives
          pNext = p+1; %Increment but leaving old in case this iteration fails
          
          augA = [zeros(m,1) A zeros(m,pNext); ones(pNext,1) zeros(pNext,n) -eye(pNext); zeros(pNext,1+n) eye(pNext)];
          augB = [b; zeros(pNext,1); repmat(UpperBoundDist,pNext,1)];         

          if isempty(lB)
              %Build bounds from maximum extents
              auglB = [0; xExtsCompact(:,1);zeros(pNext,1)];
          else
              %Build from options input
              auglB =  [0; lB;0];
          end

          if isempty(uB)
              %Build bounds from maximum extents
              auguB = [UpperBoundDist; xExtsCompact(:,2);repmat(UpperBoundDist,pNext,1)];
          else
              %Build from options input
              auguB =  [UpperBoundDist; uB;repmat(UpperBoundDist,pNext,1)];
          end
          
          %pNext
          %Prepare the distance objective function
            if pNext==1
                mComp = [];
            else
                mComp = [1 2];
            end
          NearNeighborSerialObj = @(x)-DistanceToNearestNeighbor(reshape(xMGA',1,n*pNext),pNext,n,mComp,x);
          
          problemCompact = createOptimProblem('fmincon','objective',NearNeighborSerialObj,'x0',xMGA(1,:),'Aineq',A,'bineq',b, ...
                    'Aeq',Aeq,'beq',beq,'lb',lB,'ub',uB,'nonlcon',nonlincon);
                  
          if SolveAsGA
            % Solve with a genetic algorithm
            PopSize = n*2;
            options = gaoptimset('MutationFcn',@mutationadaptfeasible,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr,@gaplotrange}, 'Display','iter','PopulationSize',PopSize,'Generations',1000,'EliteCount',ceil(0.05*PopSize),'TimeLimit',StopTime-runTime);  
            %Full formulation that incldues distance variables
               %[xNew,fBest,errorflag,Test2] = ga(@(x)-x*eye(1,1+n+pNext)',size(xMGAFull,2),augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistance(x,n,pNext,xMGA,nonlincon),options); 
            %Formulation that excludes distance variables and calculates distance in fitness function
            [xNew,fBest,errorflag,Test2] = ga(@(x)-DistanceToNearestNeighbor(x,1,n,[],xMGA),n,A,b,Aeq,beq,lB,uB,nonlincon,options); 
            sInd = 1;
            
            Test2
          else %Solve as NLP
            problemCompact.options = optimoptions('fmincon','Algorithm','sqp','maxiter',15000,'Display', 'off','MaxFunEvals',1e10);
            problemFull.options = problemCompact.options;
            
            %Try global search
            gs = GlobalSearch;
            ms = MultiStart;
            
            [xNew,fBest,errorflag,Test2] = run(gs,problemCompact);              
              
              %Calculate the new xMGAFull from prior xMGAFull alterantives.
              %Last column gets a zero!
%              xMGAFull = [xMGAFull zeros(pNext,1)];
              %xNew,fBest,errorflag,Test2] = fmincon(@(x)-x(1),sum(xMGAFull,1)/pNext,augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistance(x,n,pNext,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',100000));
               %[xNew,fBest,errorflag,Test2] = fmincon(@(x)-x*eye(1,1+n+p)',xMGAFull(1,:),augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistance(x,n,pNext,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',100000));
               %Full formulation that incldues distance variables
 %              [xNew,fBest,errorflag,Test2] = fmincon(@(x)-x(1),xMGAFull(1,:),augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistance(x,n,pNext,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',1e10));
               %Formulation that excludes distance variables and calculates
               %distance in objective function; doesn't work very well
               %[xNew,fBest,errorflag,Test2] = fmincon(@(x)-DistanceToNearestNeighbor(x,1,n,[],xMGA),xMGA(1,:),A,b,Aeq,beq,lB,uB,nonlincon,optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',150000));

              sInd = 1;
          end
          
          runTime = toc;
          
          if errorflag > 0
            %Record the new alternative
            xMGA = [xMGA;xNew(sInd:n+sInd-1)];
            %xMGAFull = [xMGAFull;xNew];
            p = pNext; %Advance the alterantive counter
            MinDist = -fBest;
            mTrialInfo = [mTrialInfo;MinDist runTime];  
              
            %Check stopping criteria
              if p>=MaxAlts
                  StopCode = 1;
              elseif runTime > StopTime
                  StopCode = 3;
              elseif MinDist <= StopDistance
                  StopCode = 2;
              end
          elseif errorflag <= -4
              StopCode = 3;
          else
              Test2.message
              StopCode = -2;             
          end
      end  
      
    elseif MGAType==3
      % Run by simultaneously searching for p maximally different
      % alternatives by maximizing the absolute value of the distance
      % between each pair of alternatives and the optimal
      % solution (L-infinity norm). Optimal solution is listed first.
       
      % Solve the problem:
      % Max Z
      % s.t Original Linear inequality, Equality, and non-linear constraints for each alt p
      %     Z <= min(sum_i(abs(x(v1,i) - x(v2,i))))   for all v1,v2, v1<v2; v2<p+1 [distance
      %                                                             measuree
      %
      % to find x(1), x(2), ..., x(p),
      % x(p+1) = the previously defined optimal solution
      %
      % note there are (p+1)*(p)/2 alternative pairs to consider in
      % calculating Z
      
      % First define the upper bound for distance variables 
      % Turn off the objfunc option before calling maxextentind/nl
      if isstruct(options) && isfield(options,'objfunc')
          optionsForExtent = rmfield(options,'objfunc');
      else
          optionsForExtent = options;
      end
      
      optionsForExtent.Algorithm = 'simplex';
      
      if blAsLP  
          [xExtsFull,xExtsCompact] = maxextentind(A,b,optionsForExtent);
          UpperBoundDist = sum(xExtsCompact(:,2)-xExtsCompact(:,1));
      else
          [xExtsFull,xExtsCompact] = maxextentindnl(A,b,optionsForExtent);
          UpperBoundDist = sum(xExtsCompact(:,2)-xExtsCompact(:,1)) ;           
      end
      
      %Formulate the augmented problem
      %Organize the decisions as (compact form) as a 1 x n*p row vector:
      %   x1(1..n) x2(i..n) .. xp(1..n)
      %Organize the constraints as:
      %  A x1  <= b      (for first alternative
      %  A x2  <= b       for second alternative
      %  :
      %  A xp  <= b      for pth alternative
      
      p = MaxAlts;
      lDistCombos = (p+1)*(p)/2;
      
      % A full version includes decision variables to represent the
      % distances between each alternative pair combo. This decision vector
      % looks like
      %  dOverall  x1(1..n) x2(i..n) .. xp(1..n)  d(xOpt,x1) d(xOpt,x2) ... d(xOpt, xp) d(x1,x2) d(x1,x3) ... d(x1,xp) ... d(xp-1,xp)
      %  
      % The full version has additional constraints for the overall distance
      % variables
      % dOverall <= d(xOpt,x1)
      % dOverall <= d(xOpt,x2)
      % :
      % dOverall <= d(xp-1,p)
      % As well as definitions of the distances
      % d(xi,xj) = sum(abs(d(xi) - d(xj))) , for all xi,xj, xj<xi, xi<=p
      % 
      
      %Build the table of paired combinations, start alternative in column
      %1, finish alternative in column 2. This table deals with the
      %lDistCombos - p combos among the alternatives themselves. The final
      %p combinations between the optimal solution and p alternatives are
      %handled later in the function that calculates the overall distance
      
      mCombos = [];
      for i=1:p-1
          mCombos = [mCombos; [repmat(i,p-i,1) [i+1:p]']];
      end
      
      [m,n] = size(A);
      
      %Original constraints for each alternative
      AeachAlt = zeros(p*m,p*n);
      for i=1:p
          AeachAlt(m*(i-1)+1:i*m,n*(i-1)+1:i*n) = A;
      end
      
      %Build the matrix inputs describing the problem in FULL form with p
      %alternatives
      augA = [zeros(p*m,1) AeachAlt zeros(p*m,lDistCombos); ones(lDistCombos,1) zeros(lDistCombos,p*n) -eye(lDistCombos); zeros(lDistCombos,1+p*n) eye(lDistCombos)];
      augB = [repmat(b,p,1); zeros(lDistCombos,1); repmat(UpperBoundDist,lDistCombos,1)];         

      if isempty(lB)
          %Build bounds from maximum extents
          auglB = [0; repmat(xExtsCompact(:,1),p,1);zeros(lDistCombos,1)];
      else
          %Build from options input
          auglB =  [0; repmat(lB,p,1);zeros(lDistCombos,1)];
      end

      if isempty(uB)
          %Build bounds from maximum extents
          auguB = [UpperBoundDist; repmat(xExtsCompact(:,2),p,1);repmat(UpperBoundDist,lDistCombos,1)];
      else
          %Build from options input
          auguB =  [UpperBoundDist; repmat(uB,p,1);repmat(UpperBoundDist,lDistCombos,1)];
      end

      %Solve for the set of alternatives in a single go.
      if IgnoreStartPoint
         xMGAStart = [];
         xMGA = zeros(1,n);
      else
         %Assume the starting point is all alterantives at the optimal value
         %(no distance)
         xMGAStart = xMGA;
      end
      xMGAFull = [0 repmat(xMGA,1,p) zeros(1,lDistCombos)];
      
      %Identify a possibly better initial starting point for the Global
      %Solver that uses the previously defined maximum extents
      %Screen out the optimal solution from the list of maximum extents
      xExtsUse = xExtsFull(any(abs(xExtsFull - repmat(xMGA,2*n,1)) > 1e-6,2),:);
      
       pToUse = min([size(xExtsUse,1) p]);
       %Reshape
       xInitToUse = reshape(xExtsUse(1:pToUse,:)',1,pToUse*n);
       
       if pToUse<p
           %Extensts were not sufficient to generate p initial alterantives. Generate remaining
           %initial alternatives to lie inside the extents. Space equally between extents along each
           %dimension
           pToAdd = p-pToUse;
           xFilled = zeros(pToAdd,n);
           pCurr = 1;
           pToFillPerDim = floor((pToAdd)/n);
           i=1;
           while (pCurr<=pToAdd) && (i<=n)
               if i==n
                   pToFillPerDim = pToAdd - pCurr + 1;
               end
               xFilled(pCurr:pCurr+pToFillPerDim-1,:)=interp1([0 1]',xExtsFull(2*i-1:2*i,:),[1:pToFillPerDim]/(1+pToFillPerDim));
               i=i+1;
               pCurr = pCurr+pToFillPerDim;
           end
           
           xInitToUse = [xInitToUse reshape(xFilled',1,n*(p-pToUse))];
       end     
       
       %Calculate the distance to the nearest neighbor as a point of
       %comparison
       dMinToNeighbor = DistanceToNearestNeighbor(xInitToUse,p,n,mCombos,xMGAStart);
       %sprintf('Minimum distance of initial solution: %d',dMinToNeighbor)
       
      %define the problem
      % Compact version
      NearNeighborObj = @(x)-DistanceToNearestNeighbor(x,p,n,mCombos,xMGAStart);      
      problemCompact = createOptimProblem('fmincon','objective',NearNeighborObj,'x0',repmat(xMGA,1,p),'Aineq',AeachAlt,'bineq',augB(1:m*p), ...
          'Aeq',Aeq,'beq',beq,'lb',auglB(2:1+n*p),'ub',auguB(2:1+n*p),'nonlcon',nonlincon);
      %Full version includes global and alternative-pair distance variables
      problemFull = createOptimProblem('fmincon','objective',@(x)-x(1),'x0',[dMinToNeighbor repmat(xMGA,1,p) repmat(dMinToNeighbor,1,lDistCombos)],'Aineq',augA,'bineq',augB, ...
          'Aeq',Aeq,'beq',beq,'lb',auglB,'ub',auguB,'nonlcon',@(x)NonLinDistanceAll(x,n,p,mCombos,xMGAStart,nonlincon));
      
      %Solve the problem
      if SolveAsGA
          % Solve with a genetic algorithm. Add required problem structure
          % elements
          PopSize = n*p;
          problemCompact.solver = 'ga';
          problemCompact.fitnessfcn = problemCompact.objective;
          problemCompact.nvars = n*p;
          problemCompact.options = gaoptimset('MutationFcn',@mutationadaptfeasible,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr,@gaplotrange}, 'Display','iter','InitialPopulation',xInitToUse,'PopulationSize',PopSize,'TimeLimit',StopTime,'Generations',1000,'EliteCount',ceil(0.05*PopSize),'Vectorized','off');  
          
          [xNew,fBest,errorflag,Test2] = ga(problemCompact);
          sInd = 1;

          %Run in vectorized form 
            %[xNew,fBest,errorflag,Test2] = ga(@(x)FittnessSimultaneous(x,p,2),size(xMGA,2),A,b,Aeq,beq,lB,uB,nonlincon,options); %
            
           Test2
       else %Solve using classic optimization. Add required problem structure elements
           problemCompact.options = optimoptions('fmincon','Algorithm','sqp','maxiter',15000,'Display', 'off','MaxFunEvals',1e10);
           problemFull.options = problemCompact.options;
           %[xNew,fNew,errorflag,Test2] = fmincon(@(x)-x*eye(1,1+p*n+lDistCombos)',xMGAFull,augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistanceAll(x,n,p,mCombos,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',100000));
           %Full form with distance variables included
            % [xNew,fBest,errorflag,Test2] = fmincon(@(x)-x(1),[dMinToNeighbor xInitToUse repmat(dMinToNeighbor,1,lDistCombos)],augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistanceAll(x,n,p,mCombos,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','interior-point','maxiter',15000,'Display', 'off','MaxFunEvals',1e10));
            %[xNew,fBest,errorflag,Test2] = fmincon(@(x)-x(1),[0 repmat(xMGA,1,p) zeros(1,lDistCombos)],augA,augB,Aeq,beq,auglB,auguB,@(x)NonLinDistanceAll(x,n,p,mCombos,xMGA,nonlincon),optimoptions(@fmincon,'Algorithm','sqp','maxiter',15000,'Display', 'off','MaxFunEvals',1e10));
            
            %Try global search
            gs = GlobalSearch;
            ms = MultiStart;
            
            [xNew,fBest,errorflag,Test2] = run(gs,problemCompact);
            dMin = DistanceToNearestNeighbor(xNew,p,n,mCombos,xMGAStart);
            %[xNew,fBest,errorflag,Test2] = fmincon(problemFull);
            sInd = 1;
       end

      if errorflag > 0
            %Record the new alternatives
            xMGA = [xMGA; reshape(xNew(sInd:p*n+sInd-1),n,p)'];
            StopCode = 1;
            MinDist = -fBest;
      elseif errorflag <=-4
            p=0;
            StopCode = 3;
      else
            Test2.message
            StopCode = -2;
            p = 0;
      end

      runTime = toc; 
    end     
    
    %Screen out the first row of xMGA which is the optimal solution
    if size(xMGA,1) > 1
        xMGA = xMGA(2:end,:);
    end
    
    StopExplain = StopCodesText{StopCode+3};
end

function [c,ceq] = NonLinDistance(x,n,p,PriorAlts,OriginalNLConstraints)
    %Generates a full set of non-linear constraints from the Original
    %non-linear contraints and distance between the current alternative x
    %and prior alternatives PriorAlts used in finding the maximum distance
    
    %INPUS
    % x = the current decision variable vector
    % n = number of regular decision variables
    % p = number of prior alternatives
    % PriorAlts = p x n matrix of the decision variable values of prior
    %       alternatives
    % OriginalConstraints = function handle to the other constraints to the
    %   original problem

    %OUTPUT
    % c = inequality constraints evaluated at x
    % ceq = equaity constraints evaluated at x
    
    ceq = [];
    c = [zeros(p,1+n) eye(p)]*x' - sum(abs(repmat(x(2:1+n),p,1)-PriorAlts),2);
    
    if ~isempty(OriginalNLConstraints)
        %embed the original constraints too
        [c_add,ceq_add] = OriginalNLConstraints(x(2:n+1));
        c = [c;c_add];
        ceq = [ceq;ceq_add];
    end
    
end

function [c,ceq] = NonLinDistanceAll(x,n,p,mOrder,xOpt,OriginalNLConstraints)
    %Generates a full set of non-linear constraints from the Original
    %non-linear contraints and distance between each pair of alternatives
    %in x and the optimal solution xOpt
    % Also, for each Original Constraint OriginalNLConstraints, generates
    % one for each alternative considered for each of the p blocks of decisions 2.. n+1, n+2.. 2n+2, ... .
    
    %INPUS
    % x = the current decision variable vector arranged as a row vector of
    %     OverallDist, Alt 1 vector 1..n, Alt 2 vector 1..n, ..., Alt p vector 1..n, Dist(1,2), Dist(1,3), ... Dist(1,p), Dist(2,3), Dist(2,4), ... Dist(2,p), ... Dist(p-1,p)
    % n = number of regular decision variables
    % p = number of alternatives
    % mOrder = lCombos x 2 matrix that specifies alternative #'s comprising the paired-combo specified by the row. First column is the start alternative. Second column is the end alternative ( 
    % xOpt = the optimal solution as a row vector (1 x n)
    % OriginalConstraints = function handle to the other constraints to the
    %   original problem

    %OUTPUT
    % c = inequality constraints evaluated at x
    % ceq = equaity constraints evaluated at x
    
    ceq = [];
    lCombos = size(mOrder,1)+p;
    % Construct matrixes that permit easily calculation of differences
    % between decision variable components among alternatives and the
    % optimal solution. First p rows compare optimal solution and the p
    % alternatives; remaining rows compare among alternatives
    xMat = reshape(x(2:p*n+1),n,p)'; %into p x n
    xFirst = [repmat(xOpt,p,1); xMat(mOrder(:,1),:)];
    xSecond = [xMat; xMat(mOrder(:,2),:)];
    
    c = [zeros(lCombos,1+n*p) eye(lCombos)]*x' - sum(abs(xFirst-xSecond),2);
    
    if ~isempty(OriginalNLConstraints)
        %embed the original constraints too, one for each block of
        %alterantives
        
        for i = 1:p
            [c_add,ceq_add] = OriginalNLConstraints(x(2+n*(i-1):n*i+1));
             c = [c;c_add];
            ceq = [ceq;ceq_add];
        end
       
    end
    
end
    

function [c,ceq] = NearOptNLConstraints(x,cFunc,rhsNearOpt,direction,OtherConstraints)
    %Generates a near-optimal constraint in the form:
    %   cFunc(x) - rhsNearOpt <= 0
    % for an optimization problem
    
    %INPUTS
    % x = the current decision variable vector
    % cFunc = a function handle to the objective function of the original
    %   underlying problem. Should be of the form cFunc(x) <= rhsNearOpt
    %   (for a minimization problem)
    % rhsNearOpt = the right hand side value of the near-optimal constraint
    %   (optimal objective function value * tolerance)
    % direction = Direction of the optimization (1=minimization,
    % -1=maximization), an adjustment applied to cFunc and rhsNearOpt
    % OtherConstraints = function handle to the other constraints to the
    %   original problem
    
    %OUTPUT
    % c = inequality constraints evaluated at x
    % ceq = equaity constraints evaluated at x
    
    ceq = [];
    c = direction*(cFunc(x) - rhsNearOpt);
    
    if nargin >=5
        %embed the original constraints too
        [c_add,ceq_add] = OtherConstraints(x);
        c = [c;c_add];
        ceq = [ceq;ceq_add];
    end
end

function [dReturn] = FittnessSimultaneous(x,p,RetPerGroup)
    %Calculates the fitness values of each alternative (row) in x
    % as the distance from the centroids of the other subpopulations
    
    % INPUT VARIABLES
    % x = q x n matrix of q alternatives (individuals) and n decision
    %    variables. x is arranged so the first 1.. q/p alternatives are
    %    part of the first sub-population, q/p+1 .. 2q/p are part of the
    %    second subpopulation, etc.
    % p = the number of desired final alternatives (sub-populations)
    % RetPerGroup = the number of alternatives to retain per
    %   sub-population. Used to rescale and re-rank fittness by
    %   sub-population
    %
    % OUTPUTS
    % dReturn = q x 1 vector describing the distance of alternative q from the
    %       centroids of each other subpopution
    %
    % Steps:
    % 1. Calculate the centroid of each subpopulation
    % 2. Calculate the distance from each alternative to the centroids of
    %    the other subpopulations
    %
    
    % Preliminaries
    [q,n] = size(x);
    if q<p
        q = p;
    else
        SubPopSize = floor(q/p);
    end
    
    SubPop = floor(([1:q]'-1)/SubPopSize)+1;
    
    if mod(q,p) > 0
        %last sub population is partial, assign to the prior subpopulation
        SubPop(end-mod(q,p):end) = p;
    end
        
    pCentroids = zeros(p,n);
    %Calculate Subpopulation Centroids
    for i=1:p
        pCentroids(i,:) = sum(x(SubPop==i,:),1)/sum(SubPop==1);
    end
    
    dInt = zeros(q,1);
    
    %Calculate distance of each individual from centroids of other
    %sub-populations
    %Create the matrixs 
    pCentroidsFull = repmat([pCentroids [1:p]'],q,1);
    pXFull = repmat([x [1:q]'],p,1);
    %Sort pXFull by the last added column so all the alternatives group
    [pXFullSort] = sortrows(pXFull,n+1);
    %Calculate distances
    dInt = sum(abs(pXFull(SubPop(pXFull(:,end))~=pCentroidsFull(:,end),1:end-1) - pCentroidsFull(SubPop(pXFull(:,end))~=pCentroidsFull(:,end),1:end-1)),2);
    dMax = max(dInt);
    
    %Find the minimum distance to another group
    dIntSh = reshape(dInt,p-1,q)';
    dMins = min(dIntSh,[],2);
    %Assign the RetPerGroup most fit members of each sub-population the same maximum
    %fitness value
    %Sort the distances by sub-population
    
    [dSort, sIndex] = sortrows([SubPop dMins],[1 2]);
    for i=1:p-1
        dSort(sIndex(SubPopSize*i + [1-RetPerGroup:0]),2) = dMax;
    end    
    dSort(end-RetPerGroup+1:end,2) = dMax;
    
    dReturn = dSort(:,2);
end

function [d]= DistanceToNearestNeighbor(x,p,n,mOrder,xOpt)
    % Computes the minimum absolute distance between a set of p alternatives
    % with n decision components and the optimal solution (-D infinity norm)
    %
    % INPUTS
    % x = 1 x n*p row vector of decision variable values for the p alternatives arrranged
    %       1..n  n+1..2n 2n+1..3n ... (p-1)n+1 .. pn 
    %       for alternatives 1,2,...,p
    % p = the number of alternatives
    % n = the number of decision variable components in each alternative
    % mOrder = an h x 2 vector of indexes representing the h paired
    %       differences between alternatives to compute. Each index
    %       references the alternative number in x
    % xOpt = 1 x n row vector of the optimal solution
    
    % reshape x into a p x n matrix where each row represents an
    % alternative
    xMat = reshape(x,n,p)';
    
    if isempty(mOrder)
        %Default distance pairs of xOpt to current
        xFirst = repmat(xMat,size(xOpt,1),1);
        xSecond = xOpt;
    else
        %Distance pairs from combinations specified in mOrder
        xFirst = xMat(mOrder(:,1),:);
        xSecond = xMat(mOrder(:,2),:);
        if ~isempty(xOpt)
            xFirst = [repmat(xOpt,p,1); xFirst];
            xSecond = [xMat; xSecond];
        end
    end
    
    d = min(sum(abs(xFirst-xSecond),2));
end