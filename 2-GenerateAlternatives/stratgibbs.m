function [X,vValid,execTime] = stratgibbs(p,problemstruct,options)
% Stratify random sample p points from within the decision space which is the convex polytope
% defined by a matrix of contraint equations Ax <= b. Stratifies to include samples
% in the polytope corners which are not ordinarily visited when uniformly
% randomly sampling high-dimensionsal polytopes.
%
% Does the strateficiation by uniformly sampling along the maximum extent
% range for each variable, then gibbs sampling within the sub-space once a
% range sample for a variable is drawn. Optionally and additionally, can
% uniform and gibbs sample along linear combinations of decision
% variables (i.e., objective functions) within the polytope (e.g., near to optimal objective function
% values) so samples span the objective space.
%
%  Method
%    1. Identify maximum extensts for each of the n variables within the
%       polytope defined by A x < b;
%    2. Uniformly randomly sample (p number of points)/(n number of variables) values from along the
%       maximum extent of the first decision variable
%    3. Take each randomly sampled value for the first decision variable,
%       reduce the polytope dimension by 1, and gibbs sample values for the remaining
%       decision variables. These values comprise one point in the
%       polytope.
%    4. Repeat Step 3 for the remaining uniformly random sampled values
%    5. Repeat Steps 3-4 for each decision variable.
%    6. If inputs include k=1:d linear combinations of decision variables, conditionally sample values for 
%          each linear combination of decision variables, then Gibbs sample
%          from the decision space using the constraints Ax <=b and c_k x =
%          f_k where f_k are the conditionally sampled values for each
%          linear combination. To determine these values, progressively build a system of d equity constraint equations
%            of the form c_k x = f_k where c_k is a vector of linear coefficients
%            for linear combination and f_k is a uniform sampled value for
%            the linear cobmination drawn from between the current maximum
%            extents of the k-th linear objective (i.e., Maxizime c_k x and
%            Minimize c_k x subject to Ax <= b, c_j x = f_j for j<k (all
%            prior linear combinations sampled)
%
%  OUTPUTS
%   X = a p-by-n matrix of random stratefied row vectors drawn
%        over the interior of the polytope
%        defined by Ax <= b. Rows 1 to p/n are points stratefied along
%        the first decision variable, rows p/n + 1 to 2*p/n are points stratified along
%        the second decision varaible, etc.
%   vValid = a p by 1 vector of integers that describe the validity of the
%         generated solution (1=valid, 0=invalid (no reason), -1=NaN,
%         -2=infeasible according to constraints)
%   execTime = execution time of the script in seconds using Matlab's tic
%         toc commands
%
%  INPUTS
%   p = number of random points to sample {single value | 2-valued vector
%         [DecisionPts LinearComboPts]). The default single value is split evenly
%         across all decision varaibles and linear combos. The dual values
%         specify the number of points among the deicsion variables and
%         linear combinations. For example, set to [0 pTs] to only sample
%         pTs points from the linear combinations.
%   problemstruct = matlab problem structure with fields for the
%   inequality, equality, lower-bound, and upper-bound constraints of the
%   problem
%      .Aineq = m x n matrix of inequality constraint coefficients
%          .Aineq x X <= .bineq
%      .bineq = m x 1 vector of right-hand side coefficients for the
%           inequality constraints (see above for .Aineq)
%      .Aeq = q x n matrix of equality constraint coeffients for the 
%           .Aeq x X = .beq
%      .beq = q x 1 vector of right hand sides of the equality constraint
%           conefficients
%      .lb = n x 1 vector of lower bounds for decision variabiles
%           X >= lb
%      .ub = n x 1 vector of upper bounds for decision variables
%           X <= ub
%
%    Also the solver .options
%
%   options
%      .matformat = {reduce [default], orig}
%          reduce = reduces A matrix rank by 1 by substracting off the
%                   sampled value before sending to the Gibbs sampler
%          orig = maintains original A matrix rank but adds two rows
%                 to represent an equity constraint for the sampled value
%                 before sending to the Gibbs sampler
%
%      .lincombo = {n by l matrix of l column vectors where each column vector contains the coefficients
%                  by which to linearly combine decision varaibles and
%                  create new projections within the polytope along which to uniformly and then
%                  Gibbs sample. (i.e., for the decision variables, these
%                  column vectors are columns from the identity matrix. For an objective function, the column
%                  represents the objective function coefficients. A multi-objective problem with nO objectives
%                  would likely have nO columns).
%
%      .errorresid = an error tolerance used to judge whether the sampled points
%               satisfy the original constraints (i.e., Ax - b <= error).
%               If a point does not, reject and resample.
%
%      .maxdrawpersample = the maximum number of attempts to Gibbs sample for a particular
%               stratefied sample. (default value = 5)
%
%      .fixeddims = a 1 x n vector of booleans that take the value of 1 to
%           indicate only stratify along the specified dimension. Default all
%           ones.
%
%      .maxtextents = n x 2 matrix of precalculated maximum extents for
%           each dimension. If omitted, the function will calculate.
%
%      .addmaxextents = takes the value of 1 to add the calculated maximum extents to
%           the alternatives returned (default value: 1 [add])
%
%      .ChainLength = the length of the Monte-Carlo Markov chain. The number of samples to draw per chain 
%           once a stratified sample is drawn and established for an objective or
%           decision variable axes. The chain represents sampling within
%           (on the interior of) the polytope. Increase to speed up
%           the overall run time as Monte Carlo Markov chain sampling is fast. (Default value: 1);
%
%      .MCMCMethod = Monte Carlo Markov Chain Method to use to sample. Will
%           call the appropriate function. Options are:
%               cprnd-gibbs {default} - a gibbs sampler written by Tim J. Benham in cprnd.m
%               cprnd-hitandrun - a hit and run sampler written by Tim J. Benham in cprnd.m
%               maxextentgibbs - a mimic of gibbs sampling that cycles through variables to find maximum extents at each step written by David E.
%                           Rosenberg in maxextentgibbs.m
% 
%      .other fields == passed along as is to the MCMC sampler
%          
% CALLED FUNCTIONS
%   - cprnd.m
%   - maxextentgibbs.m
%   - maxextentind.m
%   - linprog.m
%
%% #####################
%   Programmed by David E. Rosenberg
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%   - Febraury 2015 - Updated to specify the problem formulation as in input parameter using the matlab problem
%       structure to include inequality, inequality, lower-bound, and
%       upper-bound constraints     
%   - December 2014 - Updated to allow specifying the number of Gibbs
%       samples to draw per stratify sample (GibbsDrawsPerSample)
%   - July 2014. Updated to reject samples that are unfeasible according to the errorresid
%       criteria and add sampling capability for multiple objectives.
%   - May 8, 2013, First version
%
%   Citation:
%   David E. Rosenberg (2015) "Blended near-optimal alternative generation, 
%   visualization, and interaction for water resources decision making".
%   Water Resources Research. doi:10.1002/2013WR014667.
%   http://onlinelibrary.wiley.com/doi/10.1002/2013WR014667/full

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
    
    tStart = tic;
    n = size(a_full,2);                         % dimension
    m = size(a_full,1);                         % num constraint ineqs
    
    % Check input arguments
    
    if m < n+1                             
        % obv a prob here
        error('stratgibbs:obvprob',['Only ', num2str(m), ' inequalities. At least ',num2str(n+1),' inequalities ' ...
                            'required']);
    end
    

    %determine fixed dimensions along which to stratify sample
    if (nargin == 3) && (isstruct(options)) && isfield(options,'fixeddims')
        fixeddims = options.fixeddims;
        if size(fixeddims,2) ~= n
            warning('fixeddims is impropperly sized. Defaulting to stratify sample along all dimensions')
            fixeddims = ones(1,n);
        end
    else
        fixeddims = ones(1,n);
    end
    
    lStratDims =sum(fixeddims); %Number of dimensions along which to stratify sample
    
    ExtentAlts = [];
    AddMaxExtents = 1;
    if (nargin == 4) && (isstruct(options)) && isfield(options,'addmaxextents')
        AddMaxExtents = options.addmaxextents;
    end
    
    %Read the pre-calculated maximum extents
    PreCalcMaxExtents = [];
    if (nargin == 3) && (isstruct(options)) && isfield(options,'maxextents')
        PreCalcMaxExtents = options.maxextents;
        if size(PreCalcMaxExtents) ~= [n 2]
            warning('Precalculated maxextents are impropperly sized. Defaulting to none')
            PreCalcMaxExtents = [];
        end
        AddMaxExtents = 0;
    end
        
    %Read in the error for determining whether solutions are valid
    if (nargin == 3) && (isstruct(options)) && isfield(options,'errorresid')
        errorresid = options.errorresid;
    else
        errorresid = 0;
    end
    
    %Options structures with solve info for various routines
    if ~isfield(problemstruct,'options') || isempty(problemstruct.options)
        optLPs = struct('maxiter',1000,'Display', 'off','Algorithm','interior-point');
        problemstruct.options = optLPs;
    else
        optLPs = problemstruct.options;
    end
    
    optMEGibbs = struct('extmethod','linalg'); %opt?
    
    %Pass along the Algorithm to use
        Algorithm = optLPs.Algorithm;
        optMEGibbs.Algorithm = Algorithm;
        optLPs.Algorithm = Algorithm;
    
    %MONTE CARLO MARKOV CHAIN SAMPLER OPTIONS
    sMCMCMethods = {'cprnd-gibbs' 'cprnd-hitandrun' 'maxextentgibbs'}; 
    MCMCMethod = sMCMCMethods{1};
    
    %Read in the Monte Carlo Markov chain method
    if (nargin==3) && (isstruct(options)) && isfield(options,'MCMCMethod')
        if ismember(options.MCMCMethod,sMCMCMethods)
            MCMCMethod = options.MCMCMethod;
        else 
            warning(['stratgibbs: ', options.MCMCMethod, ' not recognized. Defaulting to ',MCMCMethod])
        end
    end
    
    ChainLength = 1;
    %Read in the Monte Carlo Markov Chain Length
    if (nargin==3) && (isstruct(options)) && isfield(options,'ChainLength')
        ChainLengthTemp = floor(options.ChainLength);
        if ChainLengthTemp < 1
            warning('GibbsDrawsPerSample must be >= 1. Proceeding with default value of 1')
        else
            ChainLength = ChainLengthTemp;
        end
    end    
    
    %determine the maxdrawpersample -- maximum number of calls to the gibbs
    %sampler per stratification sample when gibbs gives a bad sample
    if (nargin == 3) && (isstruct(options)) && isfield(options,'maxsdrawpersample')
        maxdrawpersample = options.maxdrawpersample;
    else
        maxdrawpersample = 5;
    end
    
    
    blAddCombos = 0;
    
    %determine whether there are addtional linear combos to consider
    if (nargin == 3) && (isstruct(options)) && isfield(options,'lincombo')
        lincombo = options.lincombo;
        %check the size
        [nl lLinCombos] = size(lincombo);
        
        if nl~=n %Continue with the default of no linear combinations
            warning('stratgibbs:obvprob',[num2str(n), ' decisions in A but ',num2str(nl) , ' decisions in lincombo matrix. Ignoring lincombo input']);
            pToUse = [floor(sum(p)/lStratDims)];
        else           
            %Sample a portion from the space defined by lincombos.
            %Build a new constraint matrix from the original constraints Ax <= b that define the linear
            %combinations as lLinCombos new variables in the first 1:lLinCombos columns.
            %For each linear combo add four constraints:
            %  - definition of linear combo (greater than or equal to)
            %  - definition of linear combo (less than or equal to)
            %  - global upper bound for linear combo (identified from
            %           original problem)
            %  - global lower bound (identified from the original problem)
            
            %Add the objfunc field so we can use maxextentind to identify
            %the lower and upper bounds of the linear combinations
            optionsFull.objfunc = lincombo;
            
            %[mComboFull, mComboExts] = maxextentind(A,b,optionsFull);
            [mComboFull, mComboExts] = maxextentind(problemstruct,optionsFull);
            
            lDoObjs = 1;
            
            if (any(any(isnan(mComboExts))))
                lDoObjs = 0;
            end

            %remove the lincombo field from the input structure so we avoid
            %infitite recursion
            optionsFull = rmfield(options,'lincombo');
         
            %Calculate extents again for regular decision variables
            [mExtFull, vExts] = maxextentind(problemstruct,optionsFull);
            
            if (any(any(isnan(vExts))))
                lDoObjs = 0;
            end
            
            if length(p) == 2
                pLin = round(p(2)/lLinCombos);
                pToUse = round(p(1)/lStratDims);
            else
                %Prorate by number of linear combinations and decision
                %variables
                pLin = round(sum(p)/(lLinCombos+lStratDims));
                pToUse = round(sum(p)/(lLinCombos+lStratDims));
            end
            
            %optionsFull.fixeddims = [ones(1,lLinCombos) zeros(1,n)];
            
            if isempty(PreCalcMaxExtents)
                PreCalcMaxExtentsFull = [mComboExts; vExts];
                PreCalcMaxExtents = vExts;
                if AddMaxExtents
                    ExtentAlts = mComboFull;
                end
            else
                PreCalcMaxExtentsFull = [mComboExts;PreCalcMaxExtents];
            end
            
            pMainPerLin = pLin; % floor(pLin/GibbsDrawsPerSample); ignore for objective functions
                                    
            %Call the stratify gibbs routine again with this new matrix.
            %One time for each linear combo
            Xret = []; vValidret=[];
            optionsExt=optionsFull;
            if lDoObjs==1
                
                for l=1:lLinCombos
                    
                    bRand = mComboExts(l,1) + (mComboExts(l,2)-mComboExts(l,1))*rand(pMainPerLin,1); %Random sample within maximum extents
                    CurrChainLength = 1;  % GibbsDrawsPerSample;
                    for j=1:pMainPerLin %points to randomlly sample along the stratification (a linear combination of decision variables)
                       Aeq = lincombo(:,l)';
                       Beq = bRand(j);

                       %Build up a system of equity constraint equations
                       %representing sampled linear combination values. At each
                       %step the range of allowable sampled values depends on
                       %the prior sampled value
                       if l==1
                           kInds = [2:lLinCombos];
                       elseif l==lLinCombos
                           kInds = [1:lLinCombos-1];
                       else
                           kInds = [l+1:lLinCombos 1:l-1];
                       end

                       for k=kInds
                          ProbStructCurr = problemstruct;
                          if isfield(ProbStructCurr,'Aeq')
                              ProbStructCurr.Aeq = [ProbStructCurr.Aeq;Aeq];
                              ProbStructCurr.beq = [ProbStructCurr.beq;Beq];
                          else
                              ProbStructCurr.Aeq = Aeq;
                              ProbStructCurr.beq = Beq;       
                          end
                          %Calculate the maximum extents along this dimension
                          %optionsExt.Aeq = Aeq;
                          %optionsExt.Beq = Beq;
                          optionsExt.objfunc = lincombo(:,k);

                          [xFull xCmp] = maxextentind(ProbStructCurr,optionsExt);

                          %Random sample within the extents
                          bNew = xCmp(1,1)+(xCmp(1,2)-xCmp(1,1))*rand(1,1);
                          Aeq = [Aeq;lincombo(:,k)'];
                          Beq = [Beq; bNew];                      
                       end

                       %draw a Gibbs sample from the specified linear
                       %combination values
                       if j==pMainPerLin
                           CurrChainLength = pLin - (j-1)*ChainLength;
                       end
                       
                       lDrawsNeed = CurrChainLength;
                       k=1;
                       xGood = [];
                       
                       %Update the model formulation to include the new
                       %equity constrait representing the new sampled value
                       %for the l-th linear combo (objective)
                       problemstructCurr = problemstruct;
                       if isfield(problemstructCurr,'Aeq') && isfield(problemstructCurr,'beq')
                           problemstructCurr.Aeq = [problemstructCurr.Aeq;Aeq];
                           problemstructCurr.beq = [problemstructCurr.beq;Beq];
                       else
                           problemstructCurr.Aeq = Aeq;
                           problemstructCurr.beq = Beq;
                       end
                       %Use Chebychev to draw an initial starting point
                       %(inside the inequality constraints, on the equality
                       %constriants)                      
                       x0 = chebycenterFull(problemstructCurr);
                       optMEGibbs.x0 = x0;
                       
                       if isempty(x0)
                           %Try a different formulation to generate an initial
                           %solution
                           problemstructCurr.f = circshift(eye(n,1),0)';
                           problemstructCurr.solver = 'linprog';
                           problemstructCurr.options = optLPs;

                           [x0, f_ret, f_flag] = linprog(problemstructCurr);
                           if f_flag ~= 1
                               problemstructAlt = problemstructCurr;
                               v_lb = vExts(:,1);
                               v_ub = vExts(:,2);
                               problemstructAlt.lb = v_lb;
                               problemstructAlt.ub = v_ub;
                               [xO,f_ret,f_flag] = linprog(problemstructAlt);

                               if f_flag ~= 1
                                   x0 = [];
                                   optMEGibbs.x0 = [];
                               end
                           end
                       end                     
                       
                       [a_f,b_f] = OptimiFull(problemstructCurr);
                       
                       while (lDrawsNeed>0) && (k<=maxdrawpersample)
                           %Loop to try to get good MCMC samples
 
                           %Select the Monte Carlo Markov Chain sampler
                           if strcmpi(MCMCMethod,sMCMCMethods{1})
                                Xcurr = cprnd(lDrawsNeed,a_f,b_f,struct('method','gibbs','isotropic',0,'x0',x0));                         
                           elseif strcmpi(MCMCMethod,sMCMCMethods{2})
                               Xcurr = cprnd(lDrawsNeed,a_f,b_f,struct('method','hitandrun','isotropic',0,'x0',x0));
                           elseif strcmpi(MCMCMethod,sMCMCMethods{3})
                               Xcurr = maxextentgibbs(lDrawsNeed,a_f,b_f,optMEGibbs);      
                           end
                            
                            try
                            %Not-a-number
                            vValidTemp = -1*any(isnan(Xcurr),2);
                            %Constraint feasibility
                            vResid = repmat(b_f,1,sum(vValidTemp==0)) - a_f * Xcurr';
                            vConstViolates = sum(vResid<errorresid,1)>0;
                            vValidTemp(vValidTemp==0) = -2*(vConstViolates'==1);
                            vValidTemp(vValidTemp==0) = 1;
                            catch
                               sprintf('test'); 
                            end
                            xGood = [xGood;Xcurr(vValidTemp==1,:)];
                            lDrawsNeed = CurrChainLength-size(xGood,1);
                            k=k+1;
                       end
                       
                       Xret = [Xret;Xcurr]; %store the results
                    end
                end
                
                if ~isempty(Xret)
                    lObjSols = size(Xret,1);
                    vValidret = ones(lObjSols,1);
                    %Test feasibility of Gibbs solutions
                    %for NaNs
                    vValidret(any(isnan(Xret),2)) = -1;
                    %Constraint feasibility
                    vResid = repmat(b_f,1,lObjSols) - a_f * Xret';
                    vConstViolates = sum(vResid<errorresid,1)<0;
                    vValidret(vConstViolates'==1) = -2;

                    %report feasibility statistics on these generated alternatives
                    vFeas = histc(vValidret,[-2:1]);
                    if size(vFeas,2) > 1
                        vFeas = vFeas';
                    end
                    sprintf('Feasibility of alternatives generated from linear combinations:')
                    [{'Total' 'Feasible' 'No info' 'NaNs' 'Do not satisfy constraints'}' num2cell([pLin*lLinCombos;flipud(vFeas)])]
                end
            end
            
            %Determine where to go next
            if isempty(Xret) || all(vValidret < 1)
                warning('Could not find points for linear combinations. Continuing to sample in decision space without')
                pToUse = floor(p(1)/lStratDims);
            %elseif pToUse==0
            %    %No regular sampling, we're nearly done
            %    X = Xret;
            %    vValid = vValidret;
            %    execTime = toc;
            %    return
            else
                blAddCombos = 1; %store results later            
            end
        end       
    else
        pToUse = floor(p(1)/lStratDims);
    end
    
    %We don't have to deal with linear combinations
    mLinCombos = [];
    
    %determine the matformat to use
    if (nargin == 3) && (isstruct(options)) && isfield(options,'matformat') && strcmp(lower(options.matformat),'orig')
        matformat = 'orig';
        gsize = n;
    else
        matformat = 'reduce';
        gsize = n-1;
    end
   
    %error check sample size. if only one component, it needs to be at least one
%    if (length(p)==1) && (pToUse==0)
%        warning('stratgibbs:obvprob','Need at least a sample size of one per group. Increasing the sample size to one per group.')
    if pToUse==0
        %Randomly choose a sub-set of the the decision variables along
        %which to draw a single Gibbs sample so we still get pDec samples
       vAllDIms = [1:lStratDims]';
       [uDims,bMap,cMap] = unique(cumsum(fixeddims));
         
       FixedDimsTemp = randsample(lStratDims,round(p(1)));
       fixeddimsNew = zeros(n,1);
       fixeddimsNew(bMap(FixedDimsTemp))=1; 
       fixeddims = fixeddimsNew;
       pToUse = 1;
    end
    
    %Step 1. Find independent maximum extents
    if isempty(PreCalcMaxExtents) 
        [mExtFull, vExts] = maxextentind(problemstruct,options);
        if AddMaxExtents
            ExtentAlts = [ExtentAlts;mExtFull];
        end
    else
        vExts = PreCalcMaxExtents;
    end
    
    if any(any(isnan(vExts)))
        %we had a problem
        X = []; vValid=[]; mLinCombos = [];
        execTime = [];
        return;
    end
    
    %calc the samples per decision variable stratify
    p_dec = pToUse(1);
    pMainPerDec = floor(p_dec/ChainLength);
    
    %set up the output matrix
    X = zeros((lStratDims*p_dec),n);
    vValid = zeros(lStratDims*p_dec,1);
        
    pList = [1:n]';
    %sprintf('Global Bounds on decisions\nDecision Min Max')
    %[pList vExts]
    
    mPointLogs = zeros(lStratDims,7);
    mPointLogs(1:lStratDims,1) = repmat(p_dec,lStratDims,1);
    
    %Iterate over the decision variables (Step 5)
    lStrat = 1; %count of stratify dimensions
    for i=pList(fixeddims==1)'
        %Step 2. Uniformly randomly sample 1/n-th of the points from along the
%       maximum extent of the current decision variable 
        xRndList = vExts(i,1) + (vExts(i,2)-vExts(i,1))*rand(pMainPerDec,1);
        if ChainLength==1
            X((lStrat-1)*p_dec+1:lStrat*p_dec,i) = xRndList;
        else
           %Dsitribute the sampled points into the i-th column so that for
           %GibbsDrawsPerSample>1 the values repeat on rows
           X((lStrat-1)*pMainPerDec*ChainLength+1:(lStrat)*pMainPerDec*ChainLength,i) = ... 
               reshape(repmat(xRndList',ChainLength,1),pMainPerDec*ChainLength,1);
        end
        %hist(X((i-1)*p_dec+1:i*p_dec,i));
        
        Xgibbs = zeros(ChainLength,gsize);
        c_row = circshift(eye(n,1),i-1)';
        
           for j=1:pMainPerDec    
                 
                 
               cRow = (lStrat-1)*pMainPerDec*ChainLength+(j-1)*ChainLength+1;
               cRows = [cRow:cRow+ChainLength-1];
                 
               %Update the model formulation to include the new
               %equity constrait representing the new sampled value
               %for the l-th linear combo (objective)
               problemstructCurr = problemstruct;
               if isfield(problemstructCurr,'Aeq') && isfield(problemstructCurr,'beq')
                   problemstructCurr.Aeq = [problemstructCurr.Aeq;c_row];
                   problemstructCurr.beq = [problemstructCurr.beq;X(cRow,i)];
               else
                   problemstructCurr.Aeq = c_row;
                   problemstructCurr.beq = X(cRow,i);
               end
               %Use Chebychev to draw an initial starting point
               %(inside the inequality constraints, on the equality
               %constriants)                      
               x0 = chebycenterFull(problemstructCurr);
               
               if isempty(x0)
                   %Try a different formulation to generate an initial
                   %solution
                   problemstructCurr.f = circshift(eye(n,1),i-2)';
                   problemstructCurr.solver = 'linprog';
                   problemstructCurr.options = optLPs;
                   
                   [x0, f_ret, f_flag] = linprog(problemstructCurr);
                   if f_flag ~= 1
                       problemstructAlt = problemstructCurr;
                       v_lb = vExts(:,1);
                       v_ub = vExts(:,2);
                       v_lb(i) = X(cRow,i);
                       v_ub(i) = X(cRow,i);
                       problemstructAlt.lb = v_lb;
                       problemstructAlt.ub = v_ub;
                       [xO,f_ret,f_flag] = linprog(problemstructAlt);
                       
                       if f_flag ~= 1
                           x0 = [];
                           x0red = [];
                       end
                   end
               end

               [a_f,b_f] = OptimiFull(problemstructCurr);                 

                 if strcmp('reduce',matformat)
                     %Take the m x n matrix a_f and vector b_f and reduce
                     %by the i-th column.
                     
                     a_f_red = a_f(:,i~=pList(fixeddims==1)');
                     x_curr = zeros(n,1);
                     x_curr(i) = X(cRow,i);
                     b_f_red = b_f - a_f(:,i)*x_curr(i);
                     
                     if ~isempty(x0)
                        x0red = x0(i~=pList(fixeddims==1));
                     end
                     
                     optMEGibbs.x0 = x0red;
                                          
                       %Select the Monte Carlo Markov Chain sampler
                       if strcmpi(MCMCMethod,sMCMCMethods{1})
                            Xgibbs(1:ChainLength,:) = cprnd(ChainLength,a_f_red,b_f_red,struct('method','gibbs','isotropic',0,'x0',x0red));                          
                       elseif strcmpi(MCMCMethod,sMCMCMethods{2})
                           Xgibbs(1:ChainLength,:) = cprnd(ChainLength,a_f_red,b_f_red,struct('method','hitandrun','isotropic',0,'x0',x0red));  
                       elseif strcmpi(MCMCMethod,sMCMCMethods{3})
                           Xgibbs(1:ChainLength,:) = maxextentgibbs(ChainLength,a_f_red,b_f_red,optMEGibbs);     
                       end

                     %Take the Gibbs results and put back in the larger matrix
                    if i==1
                        X(cRows,2:end) = Xgibbs(1:ChainLength,:);
                    elseif i==n
                        X(cRows,1:n-1) = Xgibbs(1:ChainLength,1:end);   
                    else
                        try
                            X(cRows,1:i-1) = Xgibbs(1:ChainLength,1:i-1);    
                            X(cRows,i+1:end) = Xgibbs(1:ChainLength,i:end);
                        catch
                            sprintf('Error')
                        end
                    end
                 else
                    if lGibbsToUse == 1
                        Xgibbs(1:ChainLength,:) = maxextentgibbs(ChainLength,a_f,b_f,optMEGibbs);
                    elseif lGibbsToUse == 2
                        Xgibbs(1:ChainLength,:) =cprnd(ChainLength,a_f,b_f,struct('method','gibbs','isotropic',0,'x0',x0));
                    end
                    X(cRows,:) = Xgibbs(1:ChainLength,:);
                 end          
            
                %Log this iteration
                mPointLogs(lStrat,4) = mPointLogs(lStrat,4)+1;
                %lGibbsDraws = lGibbsDraws+1;
                
                %Test whether stratefied/gibbs sampled solution are feasible as well as feasible according to the
                %constraints
                vNaNs = any(isnan(Xgibbs(1:ChainLength,:)),2);
                
                %Log the NaNs
                    mPointLogs(lStrat,5) = mPointLogs(lStrat,5)+sum(vNaNs); %Gibbs did not generate a solution, try again
                    vValid(cRows)=-vNaNs;
                   
                %Test feasibility of non-NaN solutions against constraints
%                 if sum(vNaNs) < lStillSample
                    dRows = cRows(vNaNs==0);
                    vResid = repmat(b_f,1,length(dRows)) - a_f * X(dRows,:)'; % Xgibbs(vNotNaNs,:);
                    vConstViolates = sum(vResid<errorresid,1)';
                    mPointLogs(i,6) = mPointLogs(i,6)+sum(vConstViolates>0); %Solution was infeasible according to the constraints, try again
                    eRows = dRows(vConstViolates > 0);
                    vValid(eRows)= -2;
                    
                %Solution is feasible, can proceed to next 
                    try
                    vValid(dRows(vConstViolates==0)) = 1;
                    catch
                        sprintf(' ')
                    end
                    mPointLogs(lStrat,7) = mPointLogs(lStrat,7)+sum(vConstViolates==0);
%                 end   
                %[i j lGibbsDraws mPointLogs(i,4) vValid((i-1)*p_dec+j)]             
           end
           
        lStrat = lStrat+1; %move to next stratum
    end
      
   if blAddCombos==1
        % Combine current results with prior linear combo results
        X = [Xret;X];
        vValid = [vValidret;vValid];
   end 
   if AddMaxExtents
        X = [ExtentAlts;X];
        vValid = [ones(size(ExtentAlts,1),1);vValid];
   end

   execTime = toc(tStart);
   
   mPointLogs(:,2) = mPointLogs(:,1) - mPointLogs(:,2);
   %mPointLogs(:,4) = mPointLogs(:,1) - mPointLogs(:,4);
   ['Target#Samples GoodInitPts OptCalls #GibbsSamples NaNGibbsSamples InfeasibleGibbsSamples FeasibleSamples']
   [pList(fixeddims==1) mPointLogs(fixeddims==1,:)]

 end
