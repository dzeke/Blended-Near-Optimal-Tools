function [X,vValid,execTime] = stratgibbs(p,A,b,options)
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
%   A = an m-by-n matrix of constraint equation
%        coefficients that includes both inequality, lower, and upper bound
%        constraints.
%   b = a m-by-1 vector of constants representing the right-hand-side of
%   the constraint equations
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
%      .Algorithm = the algorithm used to solve the underlying optimization
%           problems
%
%      .GibbsDrawsPerSample = the number of Gibbs samples to draw (interior) per
%           stratified sample drawn on a decision axis. Increase to significantly speed up
%           the overall run time as the Gibbs sampler is fast. (Default value: 1);
%
%      .other fields == passed along as is to the gibbs sampler
%          
% CALLED FUNCTIONS
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
%   - December 2014 - Updated to allow specifying the number of Gibbs
%       samples to draw per stratify sample (GibbsDrawsPerSample)
%   - July 2014. Updated to reject samples that are unfeasible according to the errorresid
%       criteria and add sampling capability for multiple objectives.
%   - May 8, 2013, First version
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

    tStart = tic;
    n = size(A,2);                         % dimension
    m = size(A,1);                         % num constraint ineqs
    
    % Check input arguments
    
    if m < n+1                             
        % obv a prob here
        error('stratgibbs:obvprob',['Only ', num2str(m), ' inequalities. At least ',num2str(n+1),' inequalities ' ...
                            'required']);
    end
    

    %determine fixed dimensions along which to stratify sample
    if (nargin == 4) && (isstruct(options)) && isfield(options,'fixeddims')
        fixeddims = options.fixeddims;
        if size(fixeddims,2) ~= n
            warning('fixeddims is impropperly sized. Defaulting to stratify sample along all dimensions')
            fixeddims = ones(1,n);
        end
    else
        fixeddims = ones(1,n);
    end
    
    lStratDims =sum(fixeddims); %Number of dimensions along which to stratify sample
    
    %Read the pre-calculated maximum extents
    PreCalcMaxExtents = [];
    if (nargin == 4) && (isstruct(options)) && isfield(options,'maxextents')
        PreCalcMaxExtents = options.maxextents;
        if size(PreCalcMaxExtents) ~= [n 2]
            warning('Precalculated maxextents are impropperly sized. Defaulting to none')
            PreCalcMaxExtents = [];
        end
    end
    
    %Read in the error for determining whether solutions are valid
    if (nargin == 4) && (isstruct(options)) && isfield(options,'errorresid')
        errorresid = options.errorresid;
    else
        errorresid = 0;
    end
    
    %Options structures with solve info for various routines
    optLPs = struct('maxiter',1000,'Display', 'off');
    optMEGibbs = struct('extmethod','linalg'); %opt?
    
    %Read in the optional Algorithm to use
    if (nargin==4) && (isstruct(options)) && isfield(options,'Algorithm')
        Algorithm = options.Algorithm;
        optMEGibbs.Algorithm = Algorithm;
        optLPs.Algorithm = Algorithm;
        optionsFull.Algorithm = Algorithm;
    end
    
    GibbsDrawsPerSample = 1;
    %Read in the GibbsDrawsPerSample
    if (nargin==4) && (isstruct(options)) && isfield(options,'GibbsDrawsPerSample')
        GibbsDrawsTemp = floor(options.GibbsDrawsPerSample);
        if GibbsDrawsTemp < 1
            warning('GibbsDrawsPerSample must be >= 1. Proceeding with default value of 1')
        else
            GibbsDrawsPerSample = GibbsDrawsTemp;
        end
    end    
    
    %determine the maxdrawpersample -- maximum number of calls to the gibbs
    %sampler per stratification sample when gibbs gives a bad sample
    if (nargin == 4) && (isstruct(options)) && isfield(options,'maxsdrawpersample')
        maxdrawpersample = options.maxdrawpersample;
    else
        maxdrawpersample = 5;
    end
    
    
    blAddCombos = 0;
    
    %determine whether there are addtional linear combos to consider
    if (nargin == 4) && (isstruct(options)) && isfield(options,'lincombo')
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
            
            [mComboFull, mComboExts] = maxextentind(A,b,optionsFull);
            
            lDoObjs = 1;
            
            if (any(any(isnan(mComboExts))))
                lDoObjs = 0;
            end

            %remove the lincombo field from the input structure so we avoid
            %infitite recursion
            optionsFull = rmfield(options,'lincombo');
         
            %Calculate extents again for regular decision variables
            [mExtFull, vExts] = maxextentind(A,b,optionsFull);
            
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
            else
                PreCalcMaxExtentsFull = [mComboExts;PreCalcMaxExtents];
            end
            
            %optionsFull.maxextents = PreCalcMaxExtentsFull;
            
            Afull = [A; lincombo';-lincombo'];
            Bfull = [b; mComboExts(:,2);-mComboExts(:,1)];
            
            pMainPerLin = pLin; % floor(pLin/GibbsDrawsPerSample); ignore for objective functions
                                    
            %Call the stratify gibbs routine again with this new matrix.
            %One time for each linear combo
            Xret = []; vValidret=[];
            optionsExt=optionsFull;
            if lDoObjs==1
                
                for l=1:lLinCombos
                    
                    bRand = mComboExts(l,1) + (mComboExts(l,2)-mComboExts(l,1))*rand(pMainPerLin,1); %Random sample within maximum extents
                    CurrGibbsDraws = 1;  % GibbsDrawsPerSample;
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
                          %Calculate the maximum extents along this dimension
                          optionsExt.Aeq = Aeq;
                          optionsExt.Beq = Beq;
                          optionsExt.objfunc = lincombo(:,k);

                          [xFull xCmp] = maxextentind(A,b,optionsExt);

                          %Random sample within the extents
                          bNew = xCmp(1,1)+(xCmp(1,2)-xCmp(1,1))*rand(1,1);
                          Aeq = [Aeq;lincombo(:,k)'];
                          Beq = [Beq; bNew];                      
                       end

                       %draw a Gibbs sample from the specified linear
                       %combination values
                       if j==pMainPerLin
                           CurrGibbsDraws = pLin - (j-1)*GibbsDrawsPerSample;
                       end
                       
                       lDrawsNeed = CurrGibbsDraws;
                       k=1;
                       xGood = [];
                       
                       while (lDrawsNeed>0) && (k<=maxdrawpersample)
                           %Loop to try to get good Gibbs samples
                            Xcurr = maxextentgibbs(lDrawsNeed,[A; Aeq; -Aeq],[b;Beq;-Beq],optMEGibbs);
                            try
                            %Not-a-number
                            vValidTemp = -1*any(isnan(Xcurr),2);
                            %Constraint feasibility
                            vResid = repmat(b,1,sum(vValidTemp==0)) - A * Xcurr';
                            vConstViolates = sum(vResid<errorresid,1)>0;
                            vValidTemp(vValidTemp==0) = -2*(vConstViolates'==1);
                            vValidTemp(vValidTemp==0) = 1;
                            catch
                               sprintf('test'); 
                            end
                            xGood = [xGood;Xcurr(vValidTemp==1,:)];
                            lDrawsNeed = CurrGibbsDraws-size(xGood,1);
                            k=k+1;
                       end
                       
                       Xret = [Xret;Xcurr]; %store the results
                    end
                end

                lObjSols = size(Xret,1);
                vValidret = ones(lObjSols,1);
                %Test feasibility of Gibbs solutions
                %for NaNs
                vValidret(any(isnan(Xret),2)) = -1;
                %Constraint feasibility
                vResid = repmat(b,1,lObjSols) - A * Xret';
                vConstViolates = sum(vResid<errorresid,1)>0;
                vValidret(vConstViolates'==1) = -2;

                %report feasibility statistics on these generated alternatives
                vFeas = histc(vValidret,[-2:1]);
                sprintf('Feasibility of alternatives generated from linear combinations:')
                [{'Total' 'Feasible' 'No info' 'NaNs' 'Do not satisfy constraints'}' num2cell([pLin*lLinCombos;flipud(vFeas)])]
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
    if (nargin == 4) && (isstruct(options)) && isfield(options,'matformat') && strcmp(lower(options.matformat),'orig')
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
        [mExtFull, vExts] = maxextentind(A,b,options);
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
    pMainPerDec = floor(p_dec/GibbsDrawsPerSample);
    
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
        if GibbsDrawsPerSample==1
            X((lStrat-1)*p_dec+1:lStrat*p_dec,i) = xRndList;
        else
           X((lStrat-1)*p_dec+1:(lStrat-1)*p_dec+pMainPerDec*GibbsDrawsPerSample,i) = ...
                reshape(repmat(xRndList',GibbsDrawsPerSample,1),pMainPerDec*GibbsDrawsPerSample,1);
           %Fill the last ones
            X((lStrat-1)*p_dec+pMainPerDec*GibbsDrawsPerSample+1:lStrat*p_dec,i) = xRndList(end);      
        end
        %hist(X((i-1)*p_dec+1:i*p_dec,i));
        
        Xgibbs = zeros(GibbsDrawsPerSample,gsize);
        c_row = circshift(eye(n,1),i-1)';
        
        %prepare the modified A matrix for the decision variable -- exclude
        %the column of the current decision variable
        if i==1
            A_mod = A(:,2:end);
        elseif i==n
            A_mod = A(:,1:end-1);
        else
            A_mod = [A(:,i+1:end) A(:,1:i-1) ]; %note, matrix is circ shifted so start on i+1 column
        end

        if strcmp('reduce',matformat)  
            %Develop the matrix to use to find a basic feasible solution to seed Gibbs
            %see if there are any all zero rows in A_mod
            v=[1:m];
            zeroRows = v(sum(A_mod==zeros(m,n-1),2)==n-1);

            if isempty(zeroRows)
                %sprintf('empty')
                A_exp = [A_mod eye(m,m-(n-1))]; %include first m-(n-1) slack variables
            else
                A_exp = A_mod;
                [mAe nAe] = size(A_exp);
                nZ = max(size(zeroRows));
                for k=1:nZ
                    A_exp = [A_exp circshift(eye(m,1),zeroRows(k)-1)];
                end
                %fill out the remaining columns to make square
                %A_exp

                if nAe+nZ < mAe
                    A_exp = [A_exp circshift(eye(m,mAe-nAe-nZ),zeroRows(nZ))];
                end              

                %A_exp
            end

            zeroRowsAft = v(sum(A_exp==zeros(m,m),2)==m);
            %size(A_exp)
            %A_exp

            %print out select variables
            if 0
               ['Cond   Rank   SparseRank']
               sprintf('%f %f %f',cond(A_exp), rank(A_exp), sprank(sparse(A_exp)))

               [A_maxRow, A_colind] = max(A_exp,[],2);
               [A_maxCol, A_rowind] = max(A_maxRow);

               LargestValue = A_exp(A_rowind,A_colind);

               A_exp(A_rowind,:) = A_exp(A_rowind,:)/LargestValue;
               sprintf('%f %f %f %f',cond(A_exp), rank(A_exp), sprank(sparse(A_exp)), LargestValue)
            end
        else
            %Prepared a second modified A matrix that includes two additional
            %rows to represent an equity constraint on the current decision
            %variable - 1st added row is the less than, 2nd added row is the
            %greater than
            %A_add = [A; c_row; -c_row];
            A_add = [circshift(A',i-1)'; eye(1,n); -eye(1,n)];
        end
        
        %xBase = circshift(eye(n,1),i-1);
               
        %A_exp_inv = inv(A_exp);
        blNeedAStart = 1;
        resid = 0;
        j=1;
        
        while (j<=p_dec) && (sum(mPointLogs(lStrat,3:4),2)<=2*maxdrawpersample*p_dec) %Step 4 (iterate over the randomly sampled values)
             %Step 3. Use the uniformly randomly sampled value for the decision variable to
             %reduce the polytope dimension by 1, and gibbs sample the remaining
             %decision variable values
             
             %xCurr = X((i-1)*p_dec+j,i)*xBase;
             %b_mod = b - A*xCurr;
             b_mod = b - A(:,i)*X((lStrat-1)*p_dec+j,i);
             
             %create the bounds associated with the equity constraint
             b_add = [b; X((lStrat-1)*p_dec+j,i); -X((lStrat-1)*p_dec+j,i)];
             
             %test whether prior starting point is feasible and we can
             %still use it to seed maxextent Gibbs
             
             if blNeedAStart == 0
                 if strcmp('reduce',matformat)
                    resid = b_mod-A_mod*x0;
                 else
                     resid = b_mod - A*x0;
                 end
             end
             
             if 1
             if (blNeedAStart) || (min(resid)<0) 
                % We need to generate a new feasible initial point
                %[mAEfin, nAEfin] = size(A_exp);
                %[X_start, fopt,errorflag] = linprog(eye(n-1,1),A_mod,b_mod,[],[],[],[],[],struct('maxiter',1000,'Display', 'off'));
                %Try optimizing. Minimize the current decision variable
                %value subject to the original constraints, current decision variable is equal to the stratified sampled value, and other decision variables within their 
                %maximum extents for other d 
                [X_start, fopt,errorflag] = linprog(circshift(eye(n,1),i),A,b,c_row,X((lStrat-1)*p_dec+j,i),vExts(1:n,1),vExts(1:n,2),[],optLPs);
             
                %Log this optimization call
                mPointLogs(lStrat,3) = mPointLogs(lStrat,3)+1;
             
                if errorflag ~= 1
                	%Couldn't find a feasibile initial point
                    Xgibbs(:,:) = NaN;
                    %['i j Min Max Sampled ErrorFlag']
                    %[i j mExtFull(2*i-1,i) mExtFull(2*i,i) X((i-1)*p_dec+j,i) errorflag]
                    
                    blNeedAStart = 1;
                    %try sampling again along the stratified decision
                    %varaible
                    X((lStrat-1)*p_dec+j,i) = vExts(i,1) + (vExts(i,2)-vExts(i,1))*rand(1,1);
                    
                    continue  %return to the top and try again
                else
                    %initial solution is feasible
                    mPointLogs(lStrat,2) = mPointLogs(lStrat,2)+1; %log as good init point
                    blNeedAStart = 0;
                    if strcmp('reduce',matformat)
                        %x0 = [X_start(1:i-1); X_start(i+1:n)];
                        x0 = [X_start(i+1:n); X_start(1:i-1)];
                    else
                        x0 = X_start;
                    end
                end
             end
             end
             
             %X_start = A_exp\b_mod; %1st row of X_start is ith+1 decision variable
             %X_start = U\(L\b_mod);
             %X_start = A_exp_inv*b_mod;
             %x0 = X_start(1:n-1);
             lGibbsDraws = 1;
             blNeedValidSol = 1;
             optMEGibbs.x0 = x0';
             
             %determine how many Gibbs samples to draw
             if j + GibbsDrawsPerSample > p_dec
                 CurrGibbsSamples = p_dec - j + 1;
             else
                 CurrGibbsSamples = GibbsDrawsPerSample;
             end
             sIndStart = (lStrat-1)*p_dec+j; %Start index
             sIndEnd = (lStrat-1)*p_dec+j + CurrGibbsSamples - 1;
             cRows = [sIndStart:sIndEnd]';
             
             while (any(vValid(cRows)<=0)) && (lGibbsDraws <= maxdrawpersample) && (j<=p_dec)  
                 
                 cStillToSample = vValid(cRows)<=0;
                 lStillSample = sum(cStillToSample);
                 
                 if strcmp('reduce',matformat)
                     %Xgibbs(j,:) = cprnd(1,A_mod,b_mod,struct('method','gibbs'));
                      %Xgibbs(j,:) = maxextentgibbs(1,A_mod,b_mod,struct('extmethod','linalg'));                     
                     Xgibbs(1:lStillSample,:) = maxextentgibbs(lStillSample,A_mod,b_mod,optMEGibbs);

                     %Take the Gibbs results and put back in the larger matrix
                    if i==1
                        X(cRows(cStillToSample),2:end) = Xgibbs(1:lStillSample,:);
                    elseif i==n
                        %size(X)
                        %size(Xgibbs)
                        %[i (i-1)*p_dec+1 i*p_dec]
                        X(cRows(cStillToSample),1:n-1) = Xgibbs(1:lStillSample,1:end);   
                    else
                        %size(X)
                        %size(Xgibbs)
                        try
                        X(cRows(cStillToSample),1:i-1) = Xgibbs(1:lStillSample,n-i+1:end);    
                        X(cRows(cStillToSample),i+1:end) = Xgibbs(1:lStillSample,1:n-i);
                        catch
                            sprintf('Error')
                        end
                    end
                 else
                    Xgibbs(1:lStillSample,:) = maxextentgibbs(lStillSample,A_add,b_add,optMEGibbs);
                    X(cRows(cStillToSample),:) = Xgibbs(1:lStillSample,:);
                 end          
            
                %Log this iteration
                mPointLogs(lStrat,4) = mPointLogs(lStrat,4)+1;
                lGibbsDraws = lGibbsDraws+1;
                
                %Test whether stratefied/gibbs sampled solution are feasible as well as feasible according to the
                %constraints
                vNaNs = any(isnan(Xgibbs(1:lStillSample,:)),2);
                
                %Log the NaNs
                    mPointLogs(lStrat,5) = mPointLogs(lStrat,5)+sum(vNaNs); %Gibbs did not generate a solution, try again
                    vValid(cRows(cStillToSample))=-vNaNs;
                   
                %Test feasibility of non-NaN solutions against constraints
                 if sum(vNaNs) < lStillSample
                    dRows = cRows(cStillToSample(vNaNs==0));
                    vResid = repmat(b,1,length(dRows)) - A * X(dRows,:)'; % Xgibbs(vNotNaNs,:);
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
                 end   
                %[i j lGibbsDraws mPointLogs(i,4) vValid((i-1)*p_dec+j)]             
             end
             
             j=j+CurrGibbsSamples; 
        end
        lStrat = lStrat+1; %move to next stratum
    end
      
   if blAddCombos==1
        % Combine current results with prior linear combo results
        X = [Xret;X];
        vValid = [vValidret;vValid];
   end   

   execTime = toc(tStart);
   
   mPointLogs(:,2) = mPointLogs(:,1) - mPointLogs(:,2);
   %mPointLogs(:,4) = mPointLogs(:,1) - mPointLogs(:,4);
   ['Target#Samples GoodInitPts OptCalls #GibbsSamples NaNGibbsSamples InfeasibleGibbsSamples FeasibleSamples']
   [pList(fixeddims==1) mPointLogs(fixeddims==1,:)]

 end
