function [X, vXs] = maxextentind(A,b,options)
% Finds and returns X -- an 2n x n matrix where each row represents a solution to the optimization problem with constraints Ax <= b. 
% Row 1 represents the lower bound on x1, row 2 the upper bound on x1, rows 3 and 4 the lower and upper bounds on x2, etc.
% Each row is independent of one another.
%
% Does this by solving 2*n optimization problems to first minimize and next
% maximize each successive decision variable
%
%  OUTPUTS
%   X = a 2n-by-n matrix where row 2i-1 is the minimum value and row 2i is the maximum value of x_i over the polytope
%        defined by AX <= b.
%   vXs = n by 2 matrix of vector bounds whose first column corresponds to the lower bound
%       and 2nd column corresponds to the upper bounds in X. i.e., row i
%       corresponds to variable i.
%
%  INPUTS
%   A = an m-by-n matrix of constraint equation
%        coefficients that includes both inequality, lower, and upper bound
%        constraints. Constraints defined in the form AX <= b.
%   b = a m-by-1 vector of constraint equation constants.
%
%   options.extmethod = method used to find the extends. There are two options
%      opt = by optimization (default)
%      linalg = by linear algebra (doesn't work)
%
%   options.objfunc = n by l matrix of l column vectors along which to
%       calculate the extents (i.e., the identity matrix would give the
%       same result as the default value) -- allows for linear combinations
%       of objectives
%
%   options.Aeq = a j by n matrix of equity constraint coefficients
%
%   options.Beq = a j by 1 vector of the right hand side values for the
%       equity constriants
%
%   options.lB = a j by 1 vector of lower bounds for the decision variables
%
%   options.uB = a j by 1 vector of upper bounds for the decision variables
%
%   options.vFixed = 1 x n vector of booleans where one means to fix the
%       the specified decision variable value at the value vFixedVal
%
%   options.vFixedVals = 1 x n vector of values to fixed the decision
%       varaible value at (when vFixed == 1). Together, vFixed and
%       vFixedVal augment rows in Aeq and Beq.
%
%   options.UnfixCurrX = boolean value 1 to indicate to unfix the constraint on the
%       current decision variable when that variable's range is being
%       indentified (keeps fixed values for other variables)
%       0 => keep all decision variable values fixed
%
%   options.Algorithm = the algorithm used to solve the underlying LP
%       problem. See linprog.m for options. Default is not specified and uses Matlab's default.
%    
%   The maximum extent method works as follows:
%    1. Start with the first decision variable (x1). Solve 2 optimization
%    problems to find both the maximum and minimum values of x1 that
%    satisfy Ax <= b. i.e.:
%
%       Max x1                         Min x1
%       such that Ax <= b       and    such that Ax < =b
%                 Aeq x = Beq                    Aeq x = Beq
%                 lB <= x <= uB                  lB <= x <= uB
%
%       to get x1(-) and x1(+)
%
%    2. Repeat for each decision variable x2, ... xn 
%
% May 1, 2013
% Updated July 2014 to include fixing decision variables at specified values and the unfixing current decision variable
% Updated Dec 2014 to allow specifying the algorithm used to solve the underlying LPs

% LICENSING
%   This code is distributed AS-IS with no expressed or implied warranty regarding the claimed functionality. The entire code or parts 
%   may be used for any non-commercial purpose so long as the use is cited per the citation above. Use for any commercial purpose requires 
%   prior written permission from the author.


% Copyright (c) 2013. David E. Rosenberg, Department of Civil and
% Environmental Engineering and Utah Water Research Lab, Utah State
% Unviersity. david.rosenberg@usu.edu

    n = size(A,2);                         % dimension
    m = size(A,1);                         % num constraint ineqs
    
    % Check input arguments
    
    if m < n+1                             
        % obv a prob here
        error('cprnd:obvprob',['Only ', num2str(m), ' inqequalities. At least ',num2str(n+1),' inequalities ' ...
                            'required']);
    end
    
    %Default values for optional paramters
    extmethod = 'opt';
    objfunc = eye(n);
    lObjs = n;
    Aeq = [];
    Beq = [];
    lB = [];
    uB = [];
    vFixed = zeros(1,n);
    vFixedVals = zeros(1,n);
    blUnfixCurrX = 0; %keep all fixed
    DfixedInd = 0; %indices of fixed variables
    AFixFull = [];
    BFixFull = [];

    %% Read Optional Inputs
    
    %determine the method to use
    if (nargin == 3) && (isstruct(options)) && isfield(options,'extmethod') && strcmp(lower(options.extmethod),'linalg')
        extmethod = 'linalg';
    end
    
    %determine if objective functions are specified
    if (nargin == 3) && (isstruct(options)) && isfield(options,'objfunc')
        objfuncTemp = options.objfunc;
        
        %verify the size of objfunc matche A
        [n2 lObjsTemp] = size(objfuncTemp);
        if n2~=n
            warning('maxextentind:obvprob',[num2str(n), 'decisions in A but ',num2str(lObjsTemp) , ' decisions in Objective function matrix. Ignoring the objective function input']);
        else
            objfunc = objfuncTemp;
            lObjs = lObjsTemp;
        end 
    end    
    
    %determine if equity constraints were passed
    if (nargin==3) && (isstruct(options)) && isfield(options,'Aeq') && isfield(options,'Beq') && ~isempty(options.Aeq) && ~isempty(options.Beq)
        AeqTemp = options.Aeq;
        BeqTemp = options.Beq;
                
        %check the sizing is correct
        if (size(AeqTemp,2) ~= n) || (size(AeqTemp,1) ~= size(BeqTemp,1))
            warning('Equity contraint parameters Aeq and Beq are improperly sized. Defaulting to no equity constraints')

        else
            Aeq=AeqTemp;
            Beq=BeqTemp;
        end
    end   
    
    %determine if lower bounds were passed
    if (nargin==3) && (isstruct(options)) && isfield(options,'lB') && ~isempty(options.lB)
        lBTemp = options.lB;
        
        %Check size is correct
        if size(lBTemp,1) ~= n
            warning('Lower bounds lB are improperly sized. Defaulting to no lower bounds')
        else
            lB = lBTemp;
        end
    end
    
    %determine if upper bounds were passed
    if (nargin==3) && (isstruct(options)) && isfield(options,'uB')  && ~isempty(options.uB)
        uBTemp = options.uB;
        
        %Check size is correct
        if size(uBTemp,1) ~= n
            warning('Upper bounds uB are improperly sized. Defaulting to no upper bounds')
        else
            uB = uBTemp;
        end
    end
    
    %Read fixed decision variable inputs
    if (nargin==3) && (isstruct(options)) && isfield(options,'vFixed') && isfield(options,'vFixedVals') && ~isempty(options.vFixed) && ~isempty(options.vFixedVals)
        vFixedTmp = options.vFixed;
        vFixedValsTmp = options.vFixedVals;
                
        %check the sizing is correct
        if (size(vFixedTmp,2) ~= n) || (size(vFixedValsTmp,2) ~= n)
            warning('vFixed and vFixedVals are improperly sized. Defaulting to no equity constraints')
        else
            %Build the Afixed and BFixed matrixes
            nDs = 1:n;
            DfixedInd = nDs(vFixedTmp==1);
            BFixFull = vFixedValsTmp(vFixedTmp==1)';
            lNumFixed = length(DfixedInd);
            for i=1:lNumFixed
                AFixFull=[AFixFull; objfunc(:,DfixedInd(i))'];
            end
            
            %Calculate the number of rows in the Aeq
            if isempty(Aeq)
                lAeqEnd  = 0;
            else
                lAeqEnd = size(Aeq,1);
            end
            
            Aeq = [Aeq; AFixFull];
            Beq = [Beq; BFixFull];
        end
    end   
    
    %determine if unfix current decision variable passed
    if (nargin==3) && (isstruct(options)) && isfield(options,'UnfixCurrX')  && ~isempty(options.UnfixCurrX)
        blUnfixCurrX = options.UnfixCurrX;
    end   
    
    %Algorithm
    if (nargin==3) && (isstruct(options)) && isfield(options,'Algorithm')
        Algorithm = options.Algorithm;
    else
        Algorithm = '';
    end
    
%% Begin Computations
    
    %extmethod;
    %set up the output extents; row 1 - min; row 2 -
    %max
    mExt = zeros(2*lObjs,n);
    vXs = zeros(lObjs,2);
    
    for j=1:lObjs
      
           %Step 1. Solve 2 optimization problems to find both the maximum
           %and minimum values of xj that satisfy Ax <= b and x(1..j-1) =
           %previously found values.
           %Create the vector of cost coefficients representing the
           %decision variable value to optimize          
           
           
           switch (extmethod)
               case 'linalg'
                   %Use linear algebra
                   Ccurr = b./(A(:,j));
                   fmin = max(Ccurr(A(:,j)<0));
                   fmax = min(Ccurr(A(:,j)>0));
                   fmaxexitflag = 1; fminexitflag = 1;
                   xmin = fmin;
                   xmax = fmax;
               case 'opt' % Use optimization
                   %f = circshift(eye(n,1),j-1);
                   f = objfunc(:,j);
                   
                   AeqToUse = Aeq;
                   BeqToUse = Beq;
                   
                   if blUnfixCurrX && ismember(j,DfixedInd)
                      %Remove the row for the current decision variable from Aeq so the decision variable limits are no longer fixed
                      RowsToUse = ~ismember([zeros(lAeqEnd,1);DfixedInd'],j);
                      AeqToUse = Aeq(RowsToUse,:);
                      BeqToUse = Beq(RowsToUse);
                   end
                   
                   LPoptions = struct('maxiter',15000,'Display', 'off');
                   if ~strcmpi(Algorithm,'')
                       LPoptions.Algorithm = Algorithm;
                   end
                   [xmin, fmin, fminexitflag minoutput] = linprog(f,A,b,AeqToUse,BeqToUse,lB,uB,[],LPoptions);
                   %if fminexitflag > 0
                      [xmax, fmax, fmaxexitflag maxoutput] = linprog(-f,A,b,AeqToUse,BeqToUse,lB,uB,[],LPoptions);
                      fmax = -fmax;
                   %end
           end
           
           if (fminexitflag>0) 
               mExt(2*j-1,:) = xmin;
               vXs(j,1) = fmin;
           else
               mExt(2*j-1,:) = NaN;
               vXs(j,1) = NaN;
           end
           
           if (fmaxexitflag>0)
              mExt(2*j,:) = xmax;  
              vXs(j,2) = fmax;
           else
              mExt(2*j,:) = NaN;
              vXs(j,2) = NaN;
           end
           
           if (fminexitflag <=0) || (fmaxexitflag <= 0)
               sprintf('j: %d, fmin: %d, fmax: %d, min flag: %d, max flag: %d',j,fmin,fmax,fminexitflag, fmaxexitflag)
               if fminexitflag <=0
                   sprintf('Min error: %s',minoutput.message)
               end
               if fmaxexitflag <=0
                   sprintf('Max error: %s',maxoutput.message)
               end
           end
    end
    
     X = mExt;
 end
