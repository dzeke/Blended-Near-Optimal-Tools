function [mValuesU nObjs vParams] = LoadEchoGamsResultsMGAComp(GDXFile,lOptSol,pIntOrig,blAddPareto,MGAMethods)
% Read in problem data and results from a GAMS GDX file, generate additional Modeling to Generate Alternatives (if specified), and prepare the results for
% plotting in the interactive near-optimal parallel coordinate plotting
% tool. Also print a table to the command window that summarizes results.
%
%   Step 1. Read in problem data and results from specified GAMS gdx file GDXFile. 
%           Results include the optimal solution, MGA alternatives generated
%           by Hop-Skip-Jump method, and Pareto solutions for a multi-objective formulation.
%           Also format everything for use in Matlab
%   Step 2. Stratify sample pIntOrig number of near-optimal alternatives using the new stratified sampling approach.
%           Calculate objective function values for the sampled alternatives
%   Step 3. Generate additional MGA alternatives in Matlab by various methods
%           (HSJ, MaxDiff-Serial, MaxDiff-Simultaneous) the user specifies in
%           the input parameter MGAMethods. 
%   Step 4. Make an interactive parallel coordinate plot for the specified groups of solutions/alternatives
%           generated in Steps 1-3. Axes including the:
%            - objective function (1st, scale at left), and 
%            - decision variables (2nd to n+1, scale at right)
%            - The Plot will have the following groups of traces
%               + Sampled near-optimal alternatives (thin light green lines)
%               + Modeling to Generate Alternatives (MGA)(thicker purple, peach,
%                   etc. lines the color corresponding to each MGA method
%                   specified in MGAMethods
%               + Optimal solution OR Pareto optimal solutions (thick black lines)
%
% INPUTS
%   - GDXFile = GAMS gdx file with optimization model outputs needed to
%       build parallel coordinate plot
%       (Optimal solution, MGA alternatives, Pareto optimal solutions, Matrix of constraint coefficients
%       and right hand sides, labels for decision variable indices -- all stored in the gams structures defined at the start of read)
%
%   - lOptSol = the row number in list of alternatives that is the optimal
%               solution to highlight on the plot with a thick black line
%
%   - pIntOrig = number of interior points to sample and add to the GAMS
%       solutions/alternatives
%
%   - blAddPareto = 1 if add additional pareto optimal alternatives at the
%       end of the alternative list and plot as a bi-objective problem.
%       The second objective is total phosphorus removed. 0 = ignore alternatives after Hop-Skip-Jump 
%       and plot as a single objective (removal cost).
% 
%   - MGAInfo = an i x 3 matrix that specifies the i method(s) to use to
%       generate MGA alternatives to include on the plot and associated
%       parameters for each method. Default value is a single scalar 0 (don't
%       include MGA methods).
%       
%          Column 1 - the MGA method to use(See further descriptions in
%                   doMGA.m)
%               0 = ignore this row
%               1 = By Hop-Skip-Jump
%               2 = Maximally different by serial
%               3 = Maximally different simultaneous (all at once)
%               4 = By Hop-Skip-Jump  run in GAMS (columns 2 and 3 ignored)
%           Column 2 - Solve by: enter 1 to solve with the Matlab genetic algorithm, 0 = with
%               classic LP or NLP optimization (Default: 0)
%           Columns 3 - maximum number of alternatives to generate
%               (Default: 10)
%
%  OUTPUTS
%   vDLabel = matrix listing of decision variable labels
%   mValuesU = m x (nObjs+nD) maxtrix of m alternatives passed to the plot function. First nObjs columns are objective function values. 
%       Next nD columns are decision variable values. 
%   nObjs = number of objectives in mValuesU
%   vParams = cell array of parameter value pairs passed to the plot function
%
%  CALLED FUNCTIONS
%    - rgdx.m (read gams gdx file)
%    - PrintGamsSo.m (to easily read in structures from GAMS gdx files)
%    - maxextentind.m (to calculate maximum extents)
%    - toolbox :  \optim\optim\linprog.m
%    - stratgibbs.m (stratify Gibbs sample from polytope)
%    - doMGA.m (generate comparison alternatives by the MGA method as
%           specified in MGAMethods)
%    - nearoptplotmo2.m (plot results in interactive parallel coordinate
%       plot tool
% 

%% #####################
%   Programmed by David E. Rosenberg
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%   - March 22, 2013. First version.
%   - Modified March 2014 to read in Pareto Solutions
%   - Modified July 2014 to read from GAMS gdx rather than numerous text outputs
%   - Modified July 2014 to read in pareto alternatives from GAMS and specify input
%       parameters for nearoptplotmo2.m to plot in multi-objectives
%   - Modified December 2014 from LoadEchoGamsResults.mga to give the user
%       full control over which Modelling to Generate Alternatives
%       method(s) is/are used as a point of comparison. See MGAMethods
%       input parameter.
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

%STEP 1 - Read in the Problem Data and Optimal Solution, MGA-HSJ, and
%Pareto Results from GAMS gdx file. Also re-format this data for use in Matlab

if ~exist(GDXFile,'file')
   error('NearOptimalLP: %s file does not exist.',GDXFile);
   return
end

%error tolerance on constraint violation calcs
epsilon = 0; %1e-4;

if nargin<5
    MGAMethods = 0;
    mgaSz=0;
else
    %Screen out rows with a method of zero   
    MGAMethods = MGAMethods(MGAMethods(:,1)>0,:);
    [mgaSz,mgaCol] = size(MGAMethods);
    
    if mgaSz ==0
        MGAMethods = 0;
        warning('No MGA methods entered. Continuing without.')
    else
        if mgaCol < 2
            %add default values for columns 2 and 3
            MGAMethods = [MGAMethods zeros(mgaSz,1) 10*ones(mgaSz,1)];
        elseif mgaCol < 3
            %add default values for column 3
            MGAMethods = [MGAMethods 10*ones(mgaSz,1)];
        end

        hMga = zeros(1,mgaSz); %Handles for the plot markers for alternatives generated by each method
    end
end

solveBy = {'Opt' 'GA'};

%Read in the scalars from the GAMS gdx file
sScalarNames = {'NumDims','GAMMA','HSJStart','NumHSJSols', 'ParetoSols' 'ExecTime'};
vScalarVals = zeros(length(sScalarNames),1);
for i = 1:length(sScalarNames);
    sScalar.name = sScalarNames{i};
    vScalarVals(i) = cell2mat(PrintGamsSo(rgdx(GDXFile,sScalar)));
end

%Read in the lengths of indexes from GAMS
sIndexNames = {'j', 'w', 's', 'i'};
vInd = zeros(length(sIndexNames),1);
cIndTEs = cell(length(sIndexNames),2);
for i=1:length(vInd)
    sName.name = sIndexNames{i};
    sName.compress = 'true';
    sName.te = 'true';
    sRes = rgdx(GDXFile,sName);
    vInd(i) = size(sRes.val,1);
    cIndTEs{i,1} = sRes.te;
    cIndTEs{i,2} = sRes.uels;
end

%Translate to variables we use in the script
nD = vScalarVals(1); %vH(3); %Number of decision variables
Gamma = vScalarVals(2); %vH(4); %Near-optimal tolerance
nJ = vInd(1); %vH(5); %Number of constraints
nW = vInd(2); %vH(6); %Number of sub-watersheds
nS = vInd(3); %vH(7); %Number of phosphorus sources
nBMP = vInd(4); %vH(8);%Number of Best Management Practices
nHSJstart = vScalarVals(3)+1; %vH(9); %Vertex number where MGA Hop-Skip-Jump alternatives start
nHSJ = vScalarVals(4); %vH(10); %Number of MGA Hop-Skip-Jump alternatives generated
nPareto = vScalarVals(5)+1; %Number of pareto alternatives genereated
GamsRunTime = vScalarVals(6); % Runtime for GAMS MGA alternatives
nParetoStart = nHSJstart+nHSJ+1; %Vertex number where pareto solutions start
nV = vScalarVals(3) + nHSJ + nPareto; %Number of vertices logged

% Include HSJ results from GAMS if MGAMethods is set to 4
blUseGamsHSJ = ismember(4,MGAMethods(:,1)); 
%               Binary takes the value of 1 to use MGA-Hop Skip Jump alternatives generated in GAMS
%                 Omitting because MGAMethods gives the user more control
%                 to generate in Matlab

%Read in Dimensioned parameters
sSingleParams = {'CAPP','ObjVal','Premove','cost','ef','AMat','vB','CTCUM'};
cSingleParams = cell(length(sSingleParams),1);
for i=1:length(sSingleParams)
    sParam.name = sSingleParams{i};
    cSingleParams{i} = PrintGamsSo(rgdx(GDXFile,sParam));   
end

%Work on CAPP to reorganize in Sub-watershed - source - BMP order (break
%out the BMPS)
cCAPPinit = cSingleParams{1};
vNum = size(cCAPPinit,1);
sNum = cell(vNum,1);
for i=1:vNum
    sTemp = strsplit(cCAPPinit{i,1},'BMP');
    sNum(i) = sTemp(2);
end
cCapp2 = [cCAPPinit num2cell(cellfun(@str2num,sNum))];
cCapp3 = sortrows(cCapp2,[2 4]);
cCapp4 = cCapp3(:,[2 1 4]);

Premove = cSingleParams{3};

%Add a column to cCapp that includes the sub-watershed
wS = unique(Premove(:,2));
cWs = [];

for i=1:nW
    cWs = [cWs; repmat(wS(i),vNum,1)];
end
cCapp5 = [cWs repmat(cCapp4,nW,1)]';

cDecHeader = cell(nD,1);
for i=1:nD
    cDecHeader{i} = strjoin(cCapp5(1:3,i)','-');
end

%Pull out the objective function values
mObjRes = cSingleParams{2};
mObjVals = mObjRes(strcmpi(mObjRes(:,2),'obj'),[1 3]);
mRemoveVals = mObjRes(strcmpi(mObjRes(:,2),'TotRem'),[1 3]);

%Transpose matrix Premove [vertex watershed source BMP PhosphorusRemoved]
%into   vertex x (Action) where action puts each Watershed-Source-BMP into
%               a separate column
PValues = zeros(nV,nD);
OValues = zeros(nV,1);
RemValues = zeros(nV,1);
vS = unique(Premove(:,1),'stable');

for i=1:size(Premove,1)
    lRow = find(strcmpi(vS,Premove(i,1)),1,'first');
    lCol = find(strcmpi(cDecHeader,strjoin(Premove(i,2:4),'-')),1,'first');
    PValues(lRow,lCol) = cell2mat(Premove(i,5));
end

for i=1:nV
    lObjRow = find(strcmpi(vS,mObjVals(i,1)),1,'first');
    OValues(lObjRow) = cell2mat(mObjVals(i,2)); 
    RemValues(lObjRow) = cell2mat(mRemoveVals(i,2));
end

%Expand the cost coefficient and efficiency values from compact to full
cVcomp = cSingleParams{4};
vEffcomp = cSingleParams{5};
cV = zeros(1,nD);
vEff = zeros(1,nD);
vLabelText = cell(nD,4);

vConstraintInds = cSingleParams{8};
lNearOptConstraint = cell2mat(vConstraintInds(4,2)); %4th is near-optimal constraint

for i=1:nD
    bmpInd = find(strcmpi(cVcomp(:,1),cCapp5(3,i)),1,'first');
    cV(i) = cell2mat(cVcomp(bmpInd,2));
    vEff(i) = cell2mat(vEffcomp(bmpInd,2));
    
    sInd = find(strcmpi(cIndTEs{3,2}{:},cCapp5(2,i)),1,'first');
    wInd = find(strcmpi(cIndTEs{2,2}{:},cCapp5(1,i)),1,'first');
    vLabelText{i,1} = cIndTEs{2,1}{wInd};
    vLabelText{i,2} = cIndTEs{3,1}{sInd};
    vLabelText{i,3} = cIndTEs{4,1}{bmpInd};
    vLabelText{i,4} = cCapp5{4,i};
end

if blAddPareto %Multi-objective
    %Combine the objective function and decision variables into a single matrix
    %with rows representing vertexes
    mValues = [OValues RemValues PValues];    
    nObjs = 2;
    vObjLabels = {'Removal Cost ($)' 'Phos. Removed (kg)'};
    
    %Add a row of all ones to the cost vector representing a second
    %objective of phosphorus removed
    cV = [cV;ones(1,nD)];
else %Single objective
    mValues = [OValues PValues]; %csvread(LPValuesFile,2,0,[2 0 1+nV nD])
    nObjs = 1;
    vObjLabels = {'Removal Cost'};
end

cVeff = (cV.*repmat(vEff,nObjs,1))';

%Read in the long descriptions of the decision variables
%[a b c d] = textread(LPDLongFile, '%s %s %s %s', 'delimiter', ',');

vLabelLong = vLabelText(:,1:3);

%vLabelLong = [a(3:end) b(3:end) c(3:end)];

%size(vLabelLong)
%substitute 'runoff' for 'difuse runoff'
for i=1:nD 
    if strcmp(vLabelLong(i,2),'Diffuse Runoff')
        vLabelLong{i,2}='Runoff';
    end
end

%Read in the contraint matrix coefficient vlaues
mA = zeros(nJ,nD);
vB = zeros(nJ,1);

vBcompact = cSingleParams{7};
for j=1:size(vBcompact,1)
    jRow = find(strcmpi(cIndTEs{1,2}{:},vBcompact(j,1)),1,'first');
    vB(jRow) = cell2mat(vBcompact(j,2));
end

mAcompact = cSingleParams{6};
for i=1:size(mAcompact,1)
    lRow = find(strcmpi(cIndTEs{1,2}{:},mAcompact(i,1)),1,'first');
    lCol = find(strcmpi(cDecHeader,strjoin(mAcompact(i,2:4),'-')),1,'first');
    mA(lRow,lCol) = cell2mat(mAcompact(i,5));   
end
%size(mTemp);


%Make a new cell array that merges all the dimensions (separated by dashes
%[-])
vLabelFull = cell(1,nD);
for i=1:nD
    %vLabelFull{i} = sprintf('%s-%s-%s',vDLabel{i,1},vDLabel{i,2},vDLabel{i,3});
    %vLabelFull{i} = sprintf('%s-%s-%s',a{i+2},b{i+2},vDLabel{i,3});
    %vLabelFull{i} = c{i+2};
    %vLabelFull{i} = sprintf('%s (%s)',c{i+2},d{i+2});
    vLabelFull{i} = sprintf('%s (B%d)',vLabelText{i,3},vLabelText{i,4});
end

%Classify each alternative into a group for later plotting
vGroupOrder = {'Near-Optimal' 'MGA-GAMS-HSJ' 'MGA-HSJ' 'MGA-Serial' 'MGA-Simultaneous' 'MGA-HSJ-Orig' 'Pareto' 'Optimum'}'; %This order ensures the optimum plots last on top and MGA on top of random interior
vLineWidth = [1 2 2 2 2 2 3 3]';

vGroupUse = zeros(1,length(vGroupOrder)); %indicator of whether group is actually used in this run
OptSolRow = [];
iGroupToHighlight = 0;

vRowCatLabel = {};
GroupSummary = []; %Summary statistics on the on the different groupings

%Optimal solutions
if lOptSol>0 
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(end),nHSJstart-1,1)];
    vGroupUse(end) = 1;
    iGroupToHighlight = vGroupOrder{end};
    OptSolRow = lOptSol;
    GroupSummary = [GroupSummary; vGroupOrder(end) solveBy(1) num2cell(nHSJstart-1) '--' '--' '--' '--' 'GAMS solve'];
end
%Hop-skip-jump alternatives from GAMS
if (nHSJ>0)
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(2),nHSJ,1)];
    if blUseGamsHSJ
        vGroupUse(2) = 1;
        GroupSummary = [GroupSummary;vGroupOrder(2) solveBy(1) num2cell(nHSJ) num2cell(GamsRunTime) num2cell(nHSJ/GamsRunTime) '--' '--' 'GAMS solve'];
        % Remove the row with the entry for MGA-HSJ from matlab
        MGAMethods = MGAMethods(MGAMethods(:,1)~=4,:);
    end
end

%Prepare plotting data to include Pareto solutions if the user specified
%this
if nPareto>0
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(end-1),nPareto,1)];
     
    if blAddPareto
        %Return and substitute pareto for the optimal group in labeling
        if lOptSol>0
            vRowCatLabel(1:nHSJstart-1) = repmat(vGroupOrder(end-1),nHSJstart-1,1);
            iGroupToHighlight = vGroupOrder{end-1};
            vGroupUse(end) = 0; %Don't use the optimal group
        end   
        
        vGroupUse(end-1) = 1;
        OptSolRow = [OptSolRow nParetoStart];
        %Reassign the line widths because so many pareto solutions
        vLineWidth = [1 1.5 1.5 1.5 1.5 1.5 2 2]';
        
        %Calculate near-optimal tolerance for second objective to allow deviations of the 2nd objective within values seen
        %when the first objective is optimized singly and separately. This is
        %simply the ratio of the 2nd objective function value obtained when the first
        %objective is optimized to the 2nd objective function value
        %obtained when the 2nd objective function is optimized.
        Gamma = [Gamma mValues(OptSolRow(1),2)/mValues(OptSolRow(2),2)];
        GroupSummary = [GroupSummary;vGroupOrder(end-1) solveBy(1) num2cell(nPareto) '--' '--' '--' '--' 'GAMS solve'];
    end
end

%Split out the objective and decision variable values
mObjs = mValues(:,1:nObjs);
mDecisions = mValues(:,nObjs+1:nObjs+nD);

%Build the matrix and vector of constraint coefficients (Ax <= b). Determine columns that always have the same value in the matrix (we will
%prune these columns to reduce dimensionality in later random sampling).
%Keep the columns whose values differ.
vPrune = [];
vKeep = [1:nD];

mAFull = mA;
vBFull = vB;
vBFullKeep = vB;

vEffFull = vEff; 
%Identify global extents for decision variables
%[mExtents, mExtCompact] = maxextentind(mAFull,vBFull);
probdef.Aineq = mAFull;
probdef.bineq = vBFull;
[mExtents,mExtCompact] = maxextentind(probdef);

%Summarize work so far
sprintf('Number of alternatives/solutions loaded from GAMS: %d',nV)
[vGroupOrder([end:-1:7 2]) num2cell([nHSJstart-1;nPareto;nHSJ])]


% STEP 2. Stratify random sample inside the polyhedron defined by the
% constraints using the new generation method.
if sum(pIntOrig) > 0
      
    %strip out the columns that all have the same values
    vMins = mExtCompact(:,1)';
    vMaxes = mExtCompact(:,2)'; %max(mDecisions);
  
    nR = max(size(vBFull));
    
    if isempty(vPrune)
        vBFullKeep = vBFull;
    else
        vBFullKeep = vBFull - mAFull(:,vPrune)*vMins(vPrune)';
    end
    
    vToUse = nV;
    mBFull = zeros(nJ,vToUse);
    for i=1:vToUse
        mBFull(:,i) = vB(:);
    end

    %Convert from decision variabls in units of phosphorus removed (kg) to area/length   
    mEffFull = repmat(vEffFull,nV,1);
    
    mAreas = mDecisions./mEffFull;
       
    %Test feasibility of the sytem of equations
    %First the near optimal vertices
    VertPts = [1:nV];
    
    mResid = (mBFull - mA*mAreas');     
    Cons = 1:nJ;
    
    %First pass, less than absolute residual    
    vConstViolates = sum(mResid<-epsilon);
    VertsThatViolate = VertPts(vConstViolates>0);
    %Second pass, less than fraction of right-hand-side value
    VertsThatViolate = VertsThatViolate(sum(mResid(:,VertsThatViolate)./mBFull(:,VertsThatViolate)<-epsilon)>0);
    
    %print out the matrix elements that are problematic
    [minResids, WorstConstraint] = min(mResid(:,VertsThatViolate));
      
    if isempty(VertsThatViolate)
        VertViolateCount = 0;
    else
        VertViolateCount = max(size(VertsThatViolate));
    end
    
    VertsToUse = VertPts(vConstViolates==0);
    %min(mResid(:,VertsThatViolate),[],1)
      
    %check that the optimal solution is the same as the first point
    
    [xopt, fopt, ferror] = linprog(cVeff(:,1),mAFull,vBFull);
    
    sprintf('Minimum objective function value %.f', min(mObjs))
    sprintf('Compare Optimal solutions\nObjective Function:')
    
    DecIndex=1:nD;
    Xdiffs = [DecIndex(xopt~=mAreas(1,:)')' xopt(DecIndex(xopt~=mAreas(1,:)')) mAreas(1,DecIndex(xopt~=mAreas(1,:)'))'];

    %Generate near-optimal alternatives by stratified gibbs sampling within
    %the near-optimal region. Exclude decision variables that are constant in the vertex set
    %Split the samples 85% / 15% between the decision variable axes and
    %objective function axes. Draw two Monte Carlo Markov chain samples per
    %group sample (chains of length 2 samples)
    
    % Full version of the model that uses decision variables of stream-bank
    % length and land area. Convert into phosphorus removed after
    % sampling. Also, puts the constraints for non-negative decision variables
    % in the Matlab .lb field which linprog handles better when
    % significantly reducing the near-optimal tolerance to 1 to generate
    % multiple optima. The original version has the non-negativity
    % constraints (x >= 0) as rows 41:79 in mAFull <= vBFullKeep (- x <= -0).
    
    %Full version
      %problStruct.Aineq = mAFull(:,vKeep);
      %problStruct.bineq = vBFullKeep;

    %Move lower-bounds into proper location
     problStruct.Aineq = mAFull(1:40,vKeep);
     vBtemp = vBFullKeep(1:40);
     %vBtemp(40) = fopt;
     problStruct.bineq = vBtemp;
     problStruct.lb = zeros(39,1);
   [NewP,vValid,sampleTime] = stratgibbs(pIntOrig*[0.85 0.15],problStruct,struct('lincombo',cVeff,'extmethod','opt', ...
            'x0',mAreas(3,vKeep)','errorresid',0,'Algorithm','interior-point','ChainLength',2,'MCMCMethod','cprnd-gibbs')); % alg = interior-point

    
    %[NewP,vValid,objvals] = stratgibbs(pInt,mAFull(:,vKeep),vBFullKeep,struct('lincombo',cVeff','extmethod','opt','x0',mAreas(3,vKeep)','errorresid',0));
   
    %eliminate invalid (infeasible) samples (shouldn't happen, but sometimes does for computational issues)
    pIntRet = size(NewP,1);   
    NewPts = [1:pIntRet];
    
    PtsThatViolate = NewPts(vValid<=0);       
    PtsToKeep = NewPts(vValid==1);
    vMatViolates = size(PtsThatViolate,2);
           
    NewP = NewP(PtsToKeep,:);
    pInt = size(NewP,1);
    
    if 0
         %Show histograms of the validity of the sampled solutions
         %Left by validity type; Right by stratify group
         fig3 = figure;
         ax1 = subplot(1,2,1,'parent',fig3);
         hist(vValid);
         hold on
         title(sprintf('vValid: %d valid',sum(vValid>0)))

         %Show which generated alternatives are in feasible
         ax2 = subplot(1,2,2,'parent',fig3);
         hist(PtsThatViolate,[(pIntRet)/(nObjs+nD)/2:pIntRet/(nObjs+nD):pIntRet])
         hold on
         title (sprintf('%d Sampled points that violate',sum(vValid<=0)));
    end

    %Add the columns back into NewP that were constant and not included in the sampling 
    NewPFull = zeros(pInt,nD);
    %NewPFull(:,vKeep) = NewP(:,:)*ef(:);
    for i=1:pInt
        NewPFull(i,vKeep) = NewP(i,:).*vEffFull(vKeep);
        NewPFull(i,vPrune) = vMins(vPrune);
    end
    
    %Count the number of sampled points that are outside the polytope
    %mResid = zeros(pInt,nD);
    mEffFull = NewPFull;    
    mBFull = zeros(nJ,pInt);

    for i=1:pInt
        mEffFull(i,:) = vEffFull(:)';
        mBFull(:,i) = vB(:);
    end
    
    lExtsToAdd = length(mExtents);
    for i=1:lExtsToAdd
      mExtents(i,:) = mExtents(i,:).*vEffFull;
    end
   
    %sprintf('Mins (Extent, Verticies, Sampled); Maxes (Extent, Verticies, Sampled)')
    %[mExtents(1,vKeep); vMins(vKeep); min(NewPFull(:,vKeep));mExtents(2,vKeep); vMaxes(vKeep); max(NewPFull(:,vKeep))]

    vDKeep = zeros(1,nD);
    vDKeep(vKeep) = vKeep;
         
    %Show some summary statistics on mResid
    VertMaxCount = zeros(2,size(vKeep,2)); %row 1=sampled pts on vertices, 2 = sampled points on global extents, 3 = verticies on global extents
    lCount = 1;
    vPts={};
    for i=1:size(vKeep,2)
        VertMaxCount(1,i) = sum(NewPFull(:,vKeep(i))>vMaxes(vKeep(i)));
        VertMaxCount(2,i) = sum(NewPFull(:,vKeep(i))>mExtents(2,vKeep(i)));
        VertMaxCount(3,i) = sum(vMaxes(vKeep(i))>mExtents(2,vKeep(i)));
        if max(VertMaxCount(3,i)) > 0
            %vPts{1,lCount}= NewPts(NewPFull(:,vKeep(i))>vMaxes(vKeep(i)));
            %vPts{3,lCount,1}= i;
            %vPts{3,lCount,2}= NewPts(vMaxes(vKeep(i))>mExtents(2,vKeep(i)));
            lCount=lCount+1;
        end
    end
    
    %Print out some statistics on the sampling
    %sprintf('Points outside:\t\t\tSampled\tVerticies\n  Matrix Constraints:\t%d\t\t%d\n  Vertex bounds:\t\t%d\t\t%d\n  Global extents:\t\t%d\t\t%d',vMatViolates,VertViolateCount,max(VertMaxCount(1,:)),0,max(VertMaxCount(2,:)),max(VertMaxCount(3,:)))
    
    %Add to the input matricies that will plot results in parallel coordinates
    %Assign each row (alternative) to a group
    vRowCatLabel = [vRowCatLabel; repmat(vGroupOrder(1),pInt,1)];
    vGroupUse(1) = 1;
    
    mDecisions = [mDecisions; NewPFull];
    %mDecisions = [mDecisions(VertsToUse,:); NewPFull];

    %Calculate objective function values for the sampled alternatives
    mObjs = [mObjs; NewPFull*cV'];
    %[nV pInt size(mObjs)]
    %mObjs = [mObjs(VertsToUse); NewPFull*cV']; 
    %sprintf('Objective Function\n  Dir.\t\tVerticies\tSampled\n  Mins\t\t%.0f\t\t%.0f\n  Maxes\t\t%.0f\t\t%.0f', ...
     %         min(mObjs(1:nV)),min(mObjs(nV+1:nV+pInt)),max(mObjs(1:nV)),max(mObjs(nV+1:nV+pInt)))
    
    %Show a histogram of the distribution of sampled objective function values
    if 0
        fig = figure;
        hist(NewPFull*cV(1,:)',[min(mObjs):(max(mObjs)-min(mObjs))/20:max(mObjs)]);
    end
    
    %Tally thes summary statistics for the group
    GroupSummaryStratSample = [vGroupOrder(1) '-' num2cell(pInt) num2cell(sampleTime) num2cell(pInt/sampleTime) num2cell(0) num2cell(0) sprintf('%d of %d feas',pIntRet-vMatViolates,pIntRet)];        
end

mValuesU = [mObjs mDecisions];

%Convert A matrix (Effective matrix) back into phosporus units nearoptplotmo2
%can use when interactively generating additional alternatives
mEff = repmat(vEffFull,nJ,1);
mAFullPhos = mAFull;
%convert but ignore the last group of non-negativity constraints;
mAFullPhos(1:vConstraintInds{4,2},:) = mAFull(1:vConstraintInds{4,2},:)./mEff(1:vConstraintInds{4,2},:);
mARet=    mAFullPhos(:,vKeep);  %mAFull(:,vKeep) ;
vBRet = vBFullKeep;

%check that the transformation to Phosphorus removed still maintains
%optimum

[xPhosOpt, fPhosOpt, fPhosError] = linprog(cV(1,:),mARet,vBRet);
fPhosOpt
[vLabelLong num2cell(xPhosOpt)]

sprintf('Minimum objective function value %.f', min(mObjs))
sprintf('Compare Optimal solutions\nObjective Function:')

%Output the Decision value and objective function ranges
sprintf('Ranges for decision variables across alternatives:')
{'Variable' 'Min' 'Max'}
[vLabelLong num2cell(mExtCompact.*[vEffFull(vKeep)' vEffFull(vKeep)'])]
sprintf('Objective Function range(s):')
num2cell([min(mObjs); max(mObjs)])

vDistLog = {};
vGroupInsert = {};

%STEP 3. GENERATE MGA Alternatives in Matlab 
% by methods (HSJ, Serial, Simultaneous) the user specifies in
% the input parameter MGAMethods. 
if mgaSz > 0
  if sum(pIntOrig) > 0 
      MGAStopTime = ceil(sampleTime);
  else
      MGAStopTime = 120;
  end
  GroupSpan = [1:max(MGAMethods(:,1))];
  GroupCount = hist(MGAMethods(:,1),GroupSpan);
  
  vLineWidthInsert = [];
  
  for i=1:size(MGAMethods,1)
    if MGAMethods(i,1) == 5
        lMGAIndex = length(vGroupOrder)-2;
    else
        lMGAIndex = 2+MGAMethods(i,1);
    end
    sprintf('Now running %s', vGroupOrder{lMGAIndex})
    %Run MGA
    if MGAMethods(i,1)==5
        % Instead run HSJ on the original model version 
        % using a different problem scaling -- e.g., decision variables in terms of
        % area/length not phosphorus removed
        [xMGAAreas, pMGA, MinDist, mIterInfo, RetFlag, RetText, rTime] = doMGA(mAreas(OptSolRow(1),:),struct('A',mAFull,'b',vBRet,'MaxAlts',MGAMethods(i,3),'StopTime',MGAStopTime,'errorcrit',1e-6,'MGAType',1,'StopDistance',0.2,'SolveAsGA',MGAMethods(i,2),'IgnoreStartPoint',0));
        %Convert from land area/stream bank length units back to phosphorus
        %removed
        if pMGA>0
            xMGA = xMGAAreas.*mEff(1:pMGA,:);        
        end
    else
        [xMGA, pMGA, MinDist, mIterInfo, RetFlag, RetText, rTime] = doMGA(mDecisions(OptSolRow(1:1+blAddPareto),:),struct('A',mARet,'b',vBRet,'MaxAlts',MGAMethods(i,3),'StopTime',MGAStopTime,'errorcrit',0,'MGAType',MGAMethods(i,1),'StopDistance',0.2,'SolveAsGA',MGAMethods(i,2),'IgnoreStartPoint',0));
    end

    %Add group labels for these alternatives
    vGroupOrderCurr = vGroupOrder(lMGAIndex);    
    
    %Add MGA alternatives to the existing alternatives
    if (RetFlag>=0) && (pMGA>0)
        mDecisions = [mDecisions; xMGA];
        %Calculate objective function values for the sampled alternatives
        mObjs = [mObjs; xMGA*cV'];
           
        %If the label has already been used, append a number to
        %differentiate and expand list of groups
        if GroupCount(MGAMethods(i,1))>1
            GroupCountTemp = hist(MGAMethods(i:end,1),GroupSpan);
            currGRep = GroupCount(MGAMethods(i,1))- GroupCountTemp(MGAMethods(i,1))+1;
            vGroupOrderCurr = {sprintf('%s%d',vGroupOrderCurr{:},currGRep)};
            vGroupInsert = [vGroupInsert;vGroupOrderCurr];
            vLineWidthInsert = [vLineWidthInsert;vLineWidth(lMGAIndex)];
        else
            vGroupUse(lMGAIndex) = 1;
        end
        vRowCatLabel = [vRowCatLabel; repmat(vGroupOrderCurr,pMGA,1)];
        %Turn the group on
          
    end

    %Update the summary statistics
    GroupSummary = [GroupSummary; vGroupOrderCurr solveBy(MGAMethods(i,2)+1) num2cell(pMGA) num2cell(rTime) num2cell(pMGA/rTime) num2cell(MinDist) num2cell(RetFlag) RetText]; 
    vDistLog{i} = mIterInfo;
  end
end

if sum(pIntOrig) > 0
    %Add the random sample group summary statistics at the end
    GroupSummary = [GroupSummary; GroupSummaryStratSample];
end

%STEP 4. PLOT UP RESULTS IN PARALLEL COORDINATES

if blAddPareto == 1 %Changes for plotting a multi-objective problem
    %Augment matrixes with one additional row representing total phosphorus
    mARet = [mARet;-ones(1,nD)];
    vBRet = [vBRet;-sum(vBRet(1:3))];
    vRowsToUse = [1:40 size(mARet,1)]';
    mAxisBounds = [900000 0 zeros(1,39); 1.1e6 16000 16000*ones(1,39)];
    lNearOptConstraint = [lNearOptConstraint; length(vRowsToUse)]';
    lShowInsetPlot = 1;
else %strip out pareto solutions and changes to plot a single-objective problem  
    mAxisBounds = [900000 zeros(1,39); 1.3e6 13000*ones(1,39)];    
    lShowInsetPlot=0;
    vRowsToUse = [1:40]';
end

%Model formulation with decision variables in units of phosphorus removed, objective function - cV
ProblemData.Aineq = mARet(vRowsToUse,:);
ProblemData.bineq =  vBRet(vRowsToUse);
ProblemData.lb = zeros(39,1);



%Insert the added MGA groups in the master list

if ~isempty(vGroupInsert)
    nGroupInsert = size(vGroupInsert,1);
    vGroupOrder = [vGroupOrder(1:end-2);vGroupInsert;vGroupOrder(end-1:end)];
    vLineWidth = [vLineWidth(1:end-2);vLineWidthInsert;vLineWidth(end-1:end)];
    vGroupUse = [vGroupUse(1:end-2) ones(1,nGroupInsert) vGroupUse(end-1:end)];
end
%Final filter of alternatives, keep alternatives associated with groups
%that have been defined in vGroupUse
vGroupShort = vGroupOrder(vGroupUse==1);
rToKeep = ismember(vRowCatLabel,vGroupShort);
rBack = cumsum(rToKeep); %Reverse mapping

%Show the results for the different MGA methods
for i=1:size(MGAMethods,1)
    if ismember(MGAMethods(i,1),[1 2 5])
        vGroupOrder{2+MGAMethods(i,1)}
        lLength = length(vDistLog{i});
        [num2cell([1:lLength]') num2cell(vDistLog{i})]
    end
end

sprintf('Runtime performance by near-optimal method')
[{'Method' 'By' 'Num. Alts.' 'Run Time (sec)' 'Rate (alts/sec)' 'Min Dist' 'Exit Code' 'Explanation'}; GroupSummary]

%Build the cell matrix of plot attributes about the groups
% Col 1 - Group Name, Col 2 - Visible, Col 3 - Line Width
mGroupData = [vGroupShort num2cell(ones(length(vGroupShort),1)) num2cell(vLineWidth(vGroupUse==1))];

%Define the optional parameters needed to plot in parallel coordinates
vParams = {'Tolerance',Gamma,'FontSize',20,'GroupToHighlight',iGroupToHighlight,'mActCat',vLabelLong(:,[1:2]),'vGroup',vRowCatLabel(rToKeep),'mGroupData',mGroupData, ...
       'vObjLabels',vObjLabels,'vXLabels',vLabelFull,'yAxisLabels',{'Removal Cost ($)' 'Phosphorus Removal (kg)'},'ProbForm',ProblemData, ... %'AMat',mARet,'Brhs', vBRet, 'cFunc',cV
       'cFunc',cV,'AxisScales','custom',[1 1],mAxisBounds,'NumTicks',5,'NearOptConstraint',lNearOptConstraint,'OptSolRow', rBack(OptSolRow)', ...
       'StartTab',1,'GenerateType',3,'GenerateMethod',2,'ShowObjsDiffColor',0,'ShowGroupLabels',1,'NumSamples',pIntOrig,'HideCheckboxes',1 ...
       'ShowControls',0,'ShowInsetPlot',lShowInsetPlot,'YAxisMargins',[0.5 0.5]};
   
%Plot with the interactive paralle coordinate tool
hNEPlot = nearoptplotmo2(mObjs(rToKeep,:), mDecisions(rToKeep,:),vParams{:});

    
    
    


