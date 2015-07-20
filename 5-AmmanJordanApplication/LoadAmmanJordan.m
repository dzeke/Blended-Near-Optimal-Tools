function LoadAmmanJordan(GDXFile,NumAlts) 
    % Read in the Amman, Jordan GAMS results from the gdx file and
    % plot in the Near-optimal parallel coordinate plotting tool
    %
    % INPUTS
    %   GDXFile = name of gdx file with GAMS results to load
    %   NumAlts = number of alternatives to read in (override number
    %       in GDX file)
    %
    %
    %  David E. Rosenberg
    %
    %  Citation:
    %  Rosenberg, D. E. (2012), Near-optimal water management to improve multi-objective decision making, 
    %  paper presented at 2012 International Congress on Environmental Modelling and Software: Managing 
    %  Resources of a Limited Planet: Pathways and Visions under Uncertainty, International Environmental 
    %  Modelling and Software Society, Leipzig, Germany. http://www.iemss.org/sites/iemss2012//proceedings/A2_0656_Rosenberg.pdf.
    
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

%%

    % Check if the GDX file exists
    if ~exist(GDXFile,'file')
       error('NearOptimalMIP: %s file does not exist.',GDXFile);
       return
    end

    %Read in the scalars from the GAMS gdx file
    sScalarNames = {'GAMMA', 'NumSols', 'NumSolvs'};
    vScalarVals = zeros(length(sScalarNames),1);
    for i = 1:length(sScalarNames);
        sScalar.name = sScalarNames{i};
        sScalarOut = rgdx(GDXFile,sScalar);
        vScalarVals(i) = sScalarOut.val;
    end

    %Translate to variables we use in the script
    Gamma = vScalarVals(1); %Near-optimal tolerance
    if nargin < 2
        NumAlts = vScalarVals(2);
    elseif NumAlts > vScalarVals(2)
        warning('Only %d alternatives in the gdx file. Ignoring the input and reading all alternatives',vScalarVals(2))
        NumAlts = vScalarVals(2);
    end
    NumSolvs = vScalarVals(3);

    %Read in the lengths of indexes from GAMS
    sIndexNames = {'i', 'cls', 'y'};
    vInd = zeros(length(sIndexNames),1);
    cIndTEs = cell(length(sIndexNames),2);
    for i=1:length(vInd)
        sName.name = sIndexNames{i};
        sName.compress = 'true';  %'true'
        sName.te = 'true';
        sRes = rgdx(GDXFile,sName);
        vInd(i) = size(sRes.val,1);
        cIndTEs{i,1} = sRes.te;
        cIndTEs{i,2} = sRes.uels;
    end

    nI = vInd(1);
    nCls = vInd(2);

    %Read in Dimensioned parameters
    sSingleParams = {'LongActs','Costs','LngCls','LevelData'};
    cSingleParams = cell(length(sSingleParams),1);
    for i=1:length(sSingleParams)
        sParam.name = sSingleParams{i};
        cSingleParams{i} = PrintGamsSo(rgdx(GDXFile,sParam));   
    end

    mCosts = cSingleParams{2};
    mActRes = cSingleParams{1};
    mLevels = cSingleParams{4};
    vActAbrev = cIndTEs{1,2}{:};

    %Reformat results so rows represent alternatives and columns attributes
    %of the alternative
    mGamsStats = cell2mat([mCosts(strcmpi(mCosts(:,4),'solstat'),6) mCosts(strcmpi(mCosts(:,4),'modstat'),6)]);
    mObjs = cell2mat([mCosts(strcmpi(mCosts(:,4),'Total'),6) mCosts(strcmpi(mCosts(:,4),'Variance'),6)]);
    mObjs = mObjs(1:NumAlts,:);
    mActTemp = cell2mat(mActRes(strcmpi(mActRes(:,5),'volum'),7));
    mActs = reshape(mActTemp(1:nI*NumAlts,:),nI,NumAlts)';   
    mLevelsOrd = reshape(cell2mat(mLevels(1:8*NumAlts,6)),8,NumAlts)';

    %Determine the validity of alternatives returned
    vModelError = (mGamsStats(:,1) == 1) + (ismember(mGamsStats(:,2),[1 2 8])) ~= 2;
    vModelError = vModelError(1:NumAlts,:);

    %Read in the objective function values
    %convert from JD to $$$
    mObjs = 1.41*mObjs; mObjs(:,2) = 1.41*mObjs(:,2);

    vActionsLong = cIndTEs{1,1};

    %Read in the action categories
    mActCls =  cSingleParams{3};
    cClassRename = [cIndTEs{2,1} {'Conservation' 'Conservation' 'New Supply' 'Conservation'}'];
    mActCats = [mActCls(:,2) mActCls(:,1)];
    for i=1:vInd(2)
        lInsts = strcmpi(mActCats(:,1),cClassRename(i,1));
        mActCats(lInsts==1,1) = repmat(cClassRename(i,2),sum(lInsts),1);
    end

    %vActionsLong= {'Meter illegal conns.','Increase meter reg.', 'Targeted cons. prog.','Drip irrigation rebates', 'Kitchen faucet rebates','Toilet rebates','Reduce physical leaks','Re-price water',  'Distant brackish water', 'Zara-Maeen project', 'Mobile desal. units', 'Red-Dead project', 'Zai plant expansion', 'Disi Conveyor', 'New local groundwater', 'New local surface water', 'Purchase Tanker Trucks',     'Reuse Wastewater' };
    %vActAbrev={'MeterIll', 'MeterReg', 'ConsProg', 'RebDIR', 'RebKFA', 'RebTDF', 'RedPhyLeak', 'RePrice', 'DesalDBW','DesalLBW','DesalMob','DesalSW','ExpandCap','NewDistGW','NewLocGW','NewSW','PurTankT','WWReuse'}

    %Create groups; identify optimal as first
    vGroupText = ['Optimal'; repmat({'Near-optimal'},NumAlts-1,1)]; 

    %Multi objective implementation
    %nearoptplotmo(0,[mData(:,1) mData(:,1)],mData(:,2:19), 1.15, zeros(1,20), zeros(1,20),ones(1,18),{'Expected cost ($ Mill/yr)' 'Cost variance ($ Mill/yr)^{2}'},vActionsLong, vActAbrev, {sprintf('%s','Objective Function Value') sprintf('%s\n%s','Implementation Level','(Million m^{3} per year)')}, 17, 1, vActCat,1)
    %single objective implementation
    %nearoptplotmo(0,mData(:,1),mData(:,2:19), 1.15, zeros(1,19), zeros(1,19),ones(1,18),{'Expected cost ($ Mill/yr)'},vActionsLong, vActAbrev, {sprintf('%s','Objective Function Value') sprintf('%s\n%s','Implementation Level','(Million m^{3} per year)')}, 22, 2, vActCat,vGroup)

    %Reorder the decision variable columns for more convienent viewing
    vActNew={'ExpandCap','ConsProg','RePrice','NewLocGW','PurTankT','DesalDBW','DesalMob','NewSW', ...   
                                                      'DesalLBW','RedPhyLeak','RebDIR','RebKFA', 'RebTDF','MeterReg','MeterIll','DesalSW','NewDistGW','WWReuse'};

    %Map the results onto the new vActNew labels provided
    uMap = MapLabels(vActAbrev,vActNew);
    mActs = mActs(:,uMap);
    vActionsLong = vActionsLong(uMap');                                                     
    mActCats = mActCats(uMap',:);
    
    %Set fixed axes and axis values for decision variables that are the
    %same throughout
    vFixed = [0 min(mActs)==max(mActs)];
    vFixedVals = [mObjs(1,1) mActs(1,:)]';
    
    %Color the axis labels in blue (New supply) and red (water
    %conservation)
    mColorBlue2Red = OSUColorRamps('BlueToRed18Step');
    mColorAxisLabels = zeros(2,1,3);
    mColorAxisLabels(1,:,:) = mColorBlue2Red(1,:);
    mColorAxisLabels(2,:,:) = flipud(mColorBlue2Red(18,:));
                                                                                                      
    %Plot with nearoptplotmo2
    vParams = {'Tolerance',Gamma,'FontSize',21,'GroupToHighlight',vGroupText{1},'mActCat',mActCats(:,1), ...
                'vFixed',vFixed,'vFixedVals',vFixedVals, ...
                'vGroup',vGroupText(vModelError==0),'vObjLabels',{'Expected Cost'}, ...
                'vXLabels',vActionsLong,'vXLabelsShort',vActNew, ...
                'yAxisLabels',{sprintf('%s\n%s','Expected Cost','($US Million)') sprintf('%s\n%s','Implementation Level','(Million m^{3} per year)')}, ...
                'AxisScales','custom',[1 1],[zeros(1,19); 50*ones(1,19)],'sGamsFile','AmmanJordanOptNearInt.gms', ...
                'GenerateType',4,'GenerateMethod',3,'mColorsAxisLabels',mColorAxisLabels};
    hNEPlot = nearoptplotmo2(mObjs(vModelError==0,1),mActs(vModelError==0,:),vParams{:});
    
    sprintf('%d alternatives loaded',NumAlts)
end