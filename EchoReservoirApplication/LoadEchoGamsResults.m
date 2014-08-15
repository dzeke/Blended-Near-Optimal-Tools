function [mValuesU nObjs vParams] = LoadEchoGamsResults(GDXFile,lOptSol,pIntOrig,blAddPareto)
% Read in GAMS results from a GDX file and prepare the results for
% plotting in the interactive near-optimal parallel coordinate plotting
% tool
%   The results include optimal solutions and maximally-different alternatives generated using Modeling to Generate Alternative (MGA)
%       hop-skip-jump methods. Additinally pareto solutions for
%       multi-objective problem.
%   The script then adds pInt # of stratified sampled near-optimal alternatives, calculates objective function values for the sampled alternatives, and 
%   makes a parallel coordinate plot with axes including the:
%   - objective function (1st, scale at left), and 
%   - decision variables (2nd to n+1, scale at right)
%   - The Plot will have the following groups of traces
%       - Random sampled interior points (thin light green lines)
%       - Modeling to Generate Alternatives (MGA) that are maximally
%         different in decision space from te optimal solution
%         (purple lines)
%       - Optimal solution OR Pareto optimal solutions (thick black lins)
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
%   - Modified July 2014 to reach in pareto alternatives and specify input
%       parameters for nearoptplotmo2.m to plot in multi-objectives
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

%STEP 1 - Read in the GAMS Optimization results

if ~exist(GDXFile,'file')
   error('NearOptimalLP: %s file does not exist.',GDXFile);
   return
end

%error tolerance on constraint violation calcs
epsilon = 0; %1e-4;

%Read in the scalars from the GAMS gdx file
sScalarNames = {'NumDims','GAMMA','HSJStart','NumHSJSols', 'ParetoSols'};
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
nV = sum(vScalarVals(3:5))+1; % min(vH(1:2)); %Max number of vertices and number of vertices used
nD = vScalarVals(1); %vH(3); %Number of decision variables
Gamma = vScalarVals(2); %vH(4); %Near-optimal tolerance
nJ = vInd(1); %vH(5); %Number of constraints
nW = vInd(2); %vH(6); %Number of sub-watersheds
nS = vInd(3); %vH(7); %Number of phosphorus sources
nBMP = vInd(4); %vH(8);%Number of Best Management Practices
nHSJstart = vScalarVals(3)+1; %vH(9); %Vertex number where MGA Hop-Skip-Jump alternatives start
nHSJ = vScalarVals(4); %vH(10); %Number of MGA Hop-Skip-Jump alternatives generated
nPareto = vScalarVals(5); %Number of pareto alternatives
nParetoStart = nHSJstart+nHSJ+1; %Vertex number where pareto solutions start

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

%STEP 2 Classify each alternative into a group for later plotting
vGroupOrder = {'Random interior' 'MGA' 'Pareto' 'Optimum'}'; %This order ensures the optimum plots last on top and MGA on top of random interior
vLineWidth = [1 2 3 3]';

vGroupUse = zeros(1,length(vGroupOrder)); %indicator of whether group is actually used in this run
OptSolRow = [];
iGroupToHighlight = 0;

vRowCatLabel = {};
%Optimal solutions
if lOptSol>0 
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(4),nHSJstart-1,1)];
    vGroupUse(4) = 1;
    iGroupToHighlight = vGroupOrder{end};
    OptSolRow = lOptSol;
end
%Hop-skip-jump alternatives
if nHSJ>0
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(2),nHSJ,1)];
    vGroupUse(2) = 1;
end
%Pareto alternatives
if nPareto>0
    vRowCatLabel = [vRowCatLabel;repmat(vGroupOrder(3),nPareto,1)];
     
    if blAddPareto
        %Return and substitute pareto for the optimal group in labeling
        if lOptSol>0
            vRowCatLabel(1:nHSJstart-1) = repmat(vGroupOrder(3),nHSJstart-1,1);
            iGroupToHighlight = vGroupOrder{3};
            vGroupUse(end) = 0; %Don't use the optimal group
        end   
        
        vGroupUse(3) = 1;
        OptSolRow = [OptSolRow; nParetoStart];
        %Reassign the line widths because so many pareto solutions
        vLineWidth = [1 1.5 2 2]';
        
        %Calculate near-optimal tolerance for second objective to allow deviations of the 2nd objective within values seen
        %when the first objective is optimized singly and separately. This is
        %simply the ratio of the 2nd objective function value obtained when the first
        %objective is optimized to the 2nd objective function value
        %obtained when the 2nd objective function is optimized.
        Gamma = [Gamma mValues(OptSolRow(1),2)/mValues(OptSolRow(2),2)];
    end
end

sprintf('Original vertices: %d',nV)
%sprintf('Original vertices: %d, Unique vertices: %d\nMGA Hop-Skip-Jump Alternatives: %d, Pareto Alternatives: %d',nV,nVU,nHSJ,nPareto)
[vGroupOrder(2:end) num2cell([nHSJ;nPareto;nHSJstart-1])]
%size(icU)
%icU

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
[mExtents, mExtCompact] = maxextentind(mAFull,vBFull);

% Stratify random sample inside the polyhedron defined by the constraints 
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
    
  
    %sum(cVeff.*xopt')
    
    %[fopt mObjs(OptVerts)']
    %[xopt mAreas(OptVerts,:)'] % mDecisions(icU(1),:)']
    DecIndex=1:nD;
    Xdiffs = [DecIndex(xopt~=mAreas(1,:)')' xopt(DecIndex(xopt~=mAreas(1,:)')) mAreas(1,DecIndex(xopt~=mAreas(1,:)'))'];

    %Generate near-optimal alternatives by stratified gibbs sampling within
    %the near-optimal region. Exclude decision variables that are constant in the vertex set
    %Split the samples 75% / 25% between the decision variable axes and
    %objective function axes
    [NewP,vValid] = stratgibbs(pIntOrig*[0.75 0.25],mAFull(:,vKeep),vBFullKeep,struct('lincombo',cVeff,'extmethod','opt','x0',mAreas(3,vKeep)','errorresid',0));
    %[NewP,vValid,objvals] = stratgibbs(pInt,mAFull(:,vKeep),vBFullKeep,struct('lincombo',cVeff','extmethod','opt','x0',mAreas(3,vKeep)','errorresid',0));
   
    %eliminate invalid samples (shouldn't happen, but sometimes does)
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
    if 1
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
    end
    
    sprintf('Points outside:\t\t\tSampled\tVerticies\n  Matrix Constraints:\t%d\t\t%d\n  Vertex bounds:\t\t%d\t\t%d\n  Global extents:\t\t%d\t\t%d',vMatViolates,VertViolateCount,max(VertMaxCount(1,:)),0,max(VertMaxCount(2,:)),max(VertMaxCount(3,:)))
    
    %Build the input matricies for plotting in parallel coordinates
    %Assign each row (alternative) to a group
    vRowCatLabel = [vRowCatLabel; repmat(vGroupOrder(1),pInt,1)];
    vGroupUse(1) = 1;
    
    mDecisions = [mDecisions; NewPFull];
    %mDecisions = [mDecisions(VertsToUse,:); NewPFull];

    %Calculate objective function values for the sampled alternatives
    mObjs = [mObjs; NewPFull*cV'];
    %[nV pInt size(mObjs)]
    %mObjs = [mObjs(VertsToUse); NewPFull*cV']; 
    sprintf('Objective Function\n  Dir.\t\tVerticies\tSampled\n  Mins\t\t%.0f\t\t%.0f\n  Maxes\t\t%.0f\t\t%.0f', ...
              min(mObjs(1:nV)),min(mObjs(nV+1:nV+pInt)),max(mObjs(1:nV)),max(mObjs(nV+1:nV+pInt)))
    
    %Show a histogram of the distribution of sampled objective function values
    if 0
        fig = figure;
        hist(NewPFull*cV(1,:)',[min(mObjs):(max(mObjs)-min(mObjs))/20:max(mObjs)]);
    end
end

%Strip out duplicate indexes
%[mValuesU, iaU, icU] = unique(mValues,'rows','stable');
%nVU = max(size(iaU));

mValuesU = [mObjs mDecisions];

%Convert A matrix back into phosporus units nearoptplotmo2
%can use when interactively generating additional alternatives
mEff = repmat(vEffFull,nJ,1);
mAFullPhos = mAFull;
%convert but ignore the non-negativity constraints;
mAFullPhos(1:40,:) = mAFull(1:40,:)./mEff(1:40,:);
mARet=    mAFullPhos(:,vKeep);  %mAFull(:,vKeep) ;
vBRet = vBFullKeep;

%Output the Decision value and objective function ranges
[vLabelLong num2cell(mExtCompact.*[vEffFull(vKeep)' vEffFull(vKeep)'])]
sprintf('Objective Function range(s):')
num2cell([min(mObjs); max(mObjs)])
%[num2cell([1:size(mObjs,1)]') vRowCatLabel]

if blAddPareto == 1 %Changes for plotting a multi-objective problem
    %Augment matrixes with one additional row representing total phosphorus
    mARet = [mARet;-ones(1,nD)];
    vBRet = [vBRet;-sum(vBRet(1:3))];
    mAxisBounds = [900000 0 zeros(1,39); 1.1e6 16000 16000*ones(1,39)];
    lNearOptConstraint = [lNearOptConstraint; length(vBRet)];
    lShowInsetPlot = 1;
else %strip out pareto solutions and changes to plot a single-objective problem  
    mAxisBounds = [900000 zeros(1,39); 1.3e6 13000*ones(1,39)];    
    lShowInsetPlot=0;
end

%Final filter of alternatives, keep alternatives associated with groups
%that have been defined in vGroupUse
vGroupShort = vGroupOrder(vGroupUse==1);
rToKeep = ismember(vRowCatLabel,vGroupShort);
rBack = cumsum(rToKeep); %Reverse mapping

%Build the cell matrix of plot attributes about the groups
% Col 1 - Group Name, Col 2 - Visible, Col 3 - Line Width
mGroupData = [vGroupShort num2cell(ones(length(vGroupShort),1)) num2cell(vLineWidth(vGroupUse==1))];

%Define the optional parameters needed to plot in parallel coordinates
vParams = {'Tolerance',Gamma,'fontsize',20,'GroupToHighlight',iGroupToHighlight,'mActCat',vLabelLong,'vGroup',vRowCatLabel(rToKeep),'mGroupData',mGroupData, ...
       'vObjLabels',vObjLabels,'vXLabels',vLabelFull,'yAxisLabels',{'Removal Cost ($)' 'Phosphorus Removal (kg)'},'AMat',mARet, ...
       'Brhs', vBRet,'cFunc',cV,'AxisScales','custom',[1 1],mAxisBounds,'NumTicks',5,'NearOptConstraint',lNearOptConstraint,'OptSolRow', rBack(OptSolRow), ...
       'StartTab',1,'GenerateType',3,'GenerateMethod',2,'ShowObjsDiffColor',0,'ShowGroupLabels',1,'NumSamples',pIntOrig,'HideCheckboxes',1 ...
       'ShowControls',0,'ShowInsetPlot',lShowInsetPlot};
   
%Plot with the interactive paralle coordinate tool
hNEPlot = nearoptplotmo2(mObjs(rToKeep,:), mDecisions(rToKeep,:),vParams{:});

    
    
    


