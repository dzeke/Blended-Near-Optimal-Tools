function [X, vXs] = maxextentind(problemstruct,options)
% Finds and returns X -- an 2n x n matrix where each row represents a solution to the optimization problem with
%   inequality, equality, lower-, and upper-bound constraints defined in
%   the problem structure problemstruct.
% Row 1 represents the lower bound on x1, row 2 the upper bound on x1, rows 3 and 4 the lower and upper bounds on x2, etc.
% Each row is independent of one another.
%
% Does this by solving 2*n optimization problems to first minimize and next
% maximize each successive decision variable. These actions represent Step
% 2 in Section 3 (Alternative Generation) in the paper Rosenberg (2015).
%
%
%  OUTPUTS
%   X = a 2n-by-n matrix where row 2i-1 is the minimum value and row 2i is the maximum value of x_i over the polytope
%        defined by AX <= b.
%   vXs = n by 2 matrix of vector bounds whose first column corresponds to the lower bound
%       and 2nd column corresponds to the upper bounds in X. i.e., row i
%       corresponds to variable i.
%
%  INPUTS
%   problemstruct = a matlab optimization problem structure with the fields
%      representing the various constraint portions of the optimization
%      problem structure (empty if omitted)
%       .Aineq = an m-by-n matrix of constraint equation
%           coefficients that includes both inequality constraints
%           defined in the form .Aineq X <= .bineq.
%
%       .bineq = a m-by-1 vector of constraint equation constants.
%
%   	.Aeq = a j by n matrix of equity constraint coefficients
%
%   	.beq = a j by 1 vector of the right hand side values for the
%           equity constriants
%
%   	.lb = a j by 1 vector of lower bounds for the decision variables
%
%   	.ub = a j by 1 vector of upper bounds for the decision variables
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
%   options.vFixed = 1 x n vector of booleans where one means to fix the
%       the specified decision variable value at the value vFixedVal
%
%   options.vFixedVals = 1 x n vector of values to fixed the decision
%       varaible value at (when vFixed == 1). Together, vFixed and
%       vFixedVal augment rows in Aeq and beq.
%
%   options.UnfixCurrX = boolean value 1 to indicate to unfix the constraint on the
%       current decision variable when that variable's range is being
%       indentified (keeps fixed values for other variables)
%       0 => keep all decision variable values fixed
%    
%   The maximum extent method works as follows:
%    1. Start with the first decision variable (x1). Solve 2 optimization
%    problems to find both the maximum and minimum values of x1 that
%    satisfy Ax <= b. i.e.:
%
%       Max x1                         Min x1
%       such that Aineq x <= bineq       and    such that Aineq x <= bineq
%                 Aeq x = beq                    Aeq x = beq
%                 lb <= x <= ub                  lb <= x <= ub
%
%       to get x1(+)                    to get x1(-)
%
%    2. Repeat for each decision variable x2, ... xn 
%
%% #####################
%   Programmed by David E. Rosenberg
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%       - May 1, 2013
%       - Updated July 2014 to include fixing decision variables at specified values and the unfixing current decision variable
%       - Updated Dec 2014 to allow specifying the algorithm used to solve the underlying LPs
%       - Updated Feb 2015 to specify problem formulation as a full Matlab structure
%
%   Citation:
%   David E. Rosenberg (2015) "Blended near-optimal alternative generation, 
%   visualization, and interaction for water resources decision making".
%   Water Resources Research. doi:10.1002/2013WR014667.
%   http://onlinelibrary.wiley.com/doi/10.1002/2013WR014667/full
%
%   LICENSING:
%   Copyright (c) 2014, David E Rosenberg
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright notice, this
%     list of conditions and the following disclaimer.
%
%   * Redistributions in binary form must reproduce the above copyright notice,
%     this list of conditions and the following disclaimer in the documentation
%     and/or other materials provided with the distribution.

%   * Neither the name of the Utah State University nor the names of its
%     contributors may be used to endorse or promote products derived from
%     this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
%   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%   Bug Reports and Feedback:
%   This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.
%
%% #######################

%   Cross from full model specification to related value
    %Check validity of full specification
    [a_full,b_full] = OptimiFull(problemstruct);
    
    if isempty(a_full)
        warning('Problem with problemstruct')
        return
    end

    Aineq = problemstruct.Aineq;
    bineq = problemstruct.bineq;

    n = size(a_full,2);%n = size(Aineq,2);                         % dimension
    m = size(a_full,1);%m = size(Aineq,1);                         % num constraint ineqs
    
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
    beq = [];
    lb = [];
    ub = [];
    vFixed = zeros(1,n);
    vFixedVals = zeros(1,n);
    blUnfixCurrX = 0; %keep all fixed
    DfixedInd = 0; %indices of fixed variables
    AFixFull = [];
    BFixFull = [];

    %% Read Optional Inputs
    
    %determine the method to use
    if (nargin == 2) && (isstruct(options)) && isfield(options,'extmethod') && strcmp(lower(options.extmethod),'linalg')
        extmethod = 'linalg';
    end
    
    %determine if objective functions are specified
    if (nargin == 2) && (isstruct(options)) && isfield(options,'objfunc')
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
    if isfield(problemstruct,'Aeq') && isfield(problemstruct,'beq') && ~isempty(problemstruct.Aeq) && ~isempty(problemstruct.beq)
        AeqTemp = problemstruct.Aeq;
        beqTemp = problemstruct.beq;
                
        %check the sizing is correct
        if (size(AeqTemp,2) ~= n) || (size(AeqTemp,1) ~= size(beqTemp,1))
            warning('Equity contraint parameters Aeq and beq are improperly sized. Defaulting to no equity constraints')

        else
            Aeq=AeqTemp;
            beq=beqTemp;
        end
    end   
    
    %determine if lower bounds were passed
    if (isstruct(problemstruct)) && isfield(problemstruct,'lb') && ~isempty(problemstruct.lb)
        lbTemp = problemstruct.lb;
        
        %Check size is correct
        if size(lbTemp,1) ~= n
            warning('Lower bounds lb are improperly sized. Defaulting to no lower bounds')
        else
            lb = lbTemp;
        end
    end
    
    %determine if upper bounds were passed
    if (isstruct(problemstruct)) && isfield(problemstruct,'ub')  && ~isempty(problemstruct.ub)
        ubTemp = problemstruct.ub;
        
        %Check size is correct
        if size(ubTemp,1) ~= n
            warning('Upper bounds ub are improperly sized. Defaulting to no upper bounds')
        else
            ub = ubTemp;
        end
    end
    
    %Read fixed decision variable inputs
    if (nargin==2) && (isstruct(options)) && isfield(options,'vFixed') && isfield(options,'vFixedVals') && ~isempty(options.vFixed) && ~isempty(options.vFixedVals)
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
            beq = [beq; BFixFull];
        end
    end   
    
    %determine if unfix current decision variable passed
    if (nargin==2) && (isstruct(options)) && isfield(options,'UnfixCurrX')  && ~isempty(options.UnfixCurrX)
        blUnfixCurrX = options.UnfixCurrX;
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
                   Ccurr = b_full./(a_full(:,j));
                   fmin = max(Ccurr(a_full(:,j)<0));
                   fmax = min(Ccurr(a_full(:,j)>0));
                   fmaxexitflag = 1; fminexitflag = 1;
                   xmin = fmin;
                   xmax = fmax;
               case 'opt' % Use optimization
                   %f = circshift(eye(n,1),j-1);
                   f = objfunc(:,j);
                   
                   AeqToUse = Aeq;
                   beqToUse = beq;
                   
                   if blUnfixCurrX && ismember(j,DfixedInd)
                      %Remove the row for the current decision variable from Aeq so the decision variable limits are no longer fixed
                      RowsToUse = ~ismember([zeros(lAeqEnd,1);DfixedInd'],j);
                      AeqToUse = Aeq(RowsToUse,:);
                      beqToUse = beq(RowsToUse);
                   end
                   
                   if isfield(problemstruct,'options') && ~isempty(problemstruct.options)
                        LPoptions = problemstruct.options;
                   else
                       LPoptions = struct('maxiter',35000,'Display','off','Algorithm','simplex');
                   end
                   [xmin, fmin, fminexitflag minoutput] = linprog(f,Aineq,bineq,AeqToUse,beqToUse,lb,ub,[],LPoptions);
                   %if fminexitflag > 0
                      [xmax, fmax, fmaxexitflag maxoutput] = linprog(-f,Aineq,bineq,AeqToUse,beqToUse,lb,ub,[],LPoptions);
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
