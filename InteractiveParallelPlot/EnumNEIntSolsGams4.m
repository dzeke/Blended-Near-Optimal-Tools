function [vObjs, mResultsInt, mResultsVal, uelsOut, vReturnFlag, mGamsStats, NumSolvs] = EnumNEIntSolsGams4(sGamsFile,vFixed,vFixedVal,uels,Tolerance,InputPassed,RunMode)
%Generate Near-Optimal Solutions to integer optimization problems
% 
% Calls and runs the GAMS file sGamsFile with decisions variables indicated by vFixed
% fixed to the values vFixedVal at the specified near optimal tolerance.
% RunMode is specific to the sGamsFile by generally specifies the type of
% solutions generated
% 
% The optimization is achieved by using the GDXMRW interface (Ferris et al, 2011, http://www.gams.com/dd/docs/tools/gdxmrw.pdf)
%  to call and execute GAMS
% 
%
% INPUT PARAMETES
% sGamsFile = string specificying the path/filename of the GAMS file to run
% vFixed = 1 x n vector of binaries  for the n decision variables in the
%           problem where a value of 1 indicates the variabled is fixed to the value
%           specified in vFixed 
% vFixedVal = 1 x n vector of values for the n decision variables in the
%           problem. The units of vFixedVal are either integer levels or
%           numbers as specified by the parameter InputPassed
% uels = 1 x n cell array of the gams uniform element labels that go with
%           vFixed and vFixedVal vectors 
% Tolerance = near optimal tolerance for the tolerable deviation constraint
%           (fraction of the optimal objective function value)
% InputPassed = integer : 1 = vFixedVal values are integer levels
%                         2 = vFixedVal values are numbers in units of the
%                        decision variable values
% RunMode = integer (value depends on settings in the sGamsFile).
%       1 = optimize only
%       2 = identify maximum extents for unfixed decision variables
%       3 = identify maximum extents for all decision variables (unfixing
%               fixed variables when their extents are evaluated)
%       4 = Enumerate all near-optimal decisions
%
% OUTPUTS
% vObjs = m x 1 vector of the objective function values for the m generated
%           near-optimal integer solutions
% mResultsInt = m x n matrix of the n decision variable levels (integer
%           values) for the m generated near-optimal integer solutions
% mResultsVal = m x n matrix of the n decision variable values (units
%            specific to the optimization program) for the m generated near-optimal
%            solutions
% uelsOut = n x 1 cell vector of the uels corresponding to the columns in
%               mResultsInt and mResultsVal. Note the ordering may be
%               differnt than the uels input depending on how GAMS returns
%               results.
% vReturnFlag = a vector of integer values that indicate the solution quality of each solution. Values
%       -2 : Error, unfeasible lower or upper bound for the provided fixed
%               inputs
%       -1 : Error, all deicison variables were fixed on input, but the
%                solution was infeasible
%        0 : Feasible lower or upper bound for the provided fixed inputs
%               but one or more decision variables still remain unfixed
%        1 : Single solution (all decision variables were fixed on input
%               and the solution was feasible
%        2 : Feasible lower or upper bound and the current decision
%               variable for which the bounds are defined was the final unfixed
%               decision variable
%        3 : Feasible lower and upper bounds on the last (unfixed) decision
%                 variable are the same so no need to further loop;
% mGamsStats = a m x 2 matrix of the gams model and solution stats for each
%           model run
% NumSolvs = the number of solves by the Gams optimization solver (both
%       that generated solutions and bounds)

% Programmed by David E. Rosenberg
% March 2014
% Utah State University
% david.rosenberg@usu.edu

    %Error checking on inputs
    if all(size(vFixed') == size(vFixedVal)) == 1
        %pivot one the inputs
        vFixed = vFixed';
    end
        
    if all(size(vFixed) == size(vFixedVal)) == 0
        error('Size of input vFixed (%d x %d) is different than vFixedVal (%d x %d)',size(vFixed), size(vFixedVal))
        return
    end
    
    if size(vFixed,2)==1
        warning('vFixed and vFixedVals are column vectors. Converting to row vectors')
        vFixed = vFixed';
        vFixedVal = vFixedVal';
    end

    %Echo the input

    %Tolerance
    %InputPassed
    
    
    %Initialize Return values
    vObjs = [];
    mResultsInt = [];
    mResultsVal = [];
    vReturnFlag = [];
    uelsOut = {};
    mGamsStats = [];
    NumSolvs = 0;
        
    n = size(vFixed,2);

    %Prepare inputs for GAMS 
    sLFixed.name = 'LFixed';
    sLFixed.val = vFixed;
    sLFixed.uels = uels;
    sLFixed.type = 'parameter';
    sLFixed.form = 'full';
    
    sInputPassed.name = 'InputPassed';
    sInputPassed.type = 'parameter';    
    sInputPassed.val = InputPassed;
    
    sLFixedVal.name = 'LFixedVal';
    sLFixedVal.val = vFixedVal;
    sLFixedVal.uels = uels;
    sLFixedVal.type = 'parameter';
    sLFixedVal.form = 'full';
    
    sLFixedValVol.name = 'LFixedValVol';    
    sLFixedValVol.val = vFixedVal;
    sLFixedValVol.uels = uels;
    sLFixedValVol.type = 'parameter';
    sLFixedValVol.form = 'full';
    
    sRunMode.name = 'RunMode';
    sRunMode.type = 'parameter';    
    sRunMode.val = RunMode;        

    gamso.input = 'exec';
    gamso.output = '';
           
        sGAMMA.name = 'GAMMA';
        sGAMMA.val = Tolerance;
        sGAMMA.type = 'parameter';
       
   %Write the input data to GDX
   wgdx('matdata2.gdx',sLFixed,sLFixedVal,sLFixedValVol,sInputPassed,sGAMMA,sRunMode)
   
   %Run gams
   gams(sGamsFile);
   
   %Read from the output GDX file
   sCosts.name = 'Costs';
   sLongActs.name = 'LongActs';
   sLevels.name = 'LevelData';
   sLToUse.name='LToUse';
   sNumSols.name='NumSols';
   sNumSolvs.name='NumSolvs';
   sTurnedFixed.name='vTurnedFixed';
   
   mCosts = PrintGamsSo(rgdx('matsol.gdx',sCosts));
   mActRes = PrintGamsSo(rgdx('matsol.gdx',sLongActs));
   mLevels = PrintGamsSo(rgdx('matsol.gdx',sLevels));
   mLToUse = PrintGamsSo(rgdx('matsol.gdx',sLToUse));
   mTurnedFixed = PrintGamsSo(rgdx('matsol.gdx',sTurnedFixed));
   sNumSols = rgdx('matsol.gdx',sNumSols);
   sNumSolves = rgdx('matsol.gdx',sNumSolvs);
   
   NumSolutions = sNumSols.val;
   NumSolvs = sNumSolves.val;
   
   uelsOut = mLToUse(:,1);

   mGamsStats = cell2mat([mCosts(strcmpi(mCosts(:,4),'solstat'),6) mCosts(strcmpi(mCosts(:,4),'modstat'),6)]);
   vObjs = cell2mat(mCosts(strcmpi(mCosts(:,4),'Total'),6));		
   mResultsInt = reshape(cell2mat(mActRes(strcmpi(mActRes(:,5),'level'),7)),n,NumSolutions)';
   mResultsVal = reshape(cell2mat(mActRes(strcmpi(mActRes(:,5),'volum'),7)),n,NumSolutions)';
   
   mLevelsOrd = reshape(cell2mat(mLevels(:,6)),7,NumSolutions)';
      
   %sprintf('The following decision variables had the same lower and upper maximal extents and were fixed in the analysis at the following values:')
   %DecResults = [mTurnedFixed mResultsVal(1,strcmpi(mTurnedFixed(:,1),uelsOut))'];
   %DecResults(cell2mat(DecResults(:,2))==1,[1 3])
			
   vReturnFlag = mLevelsOrd(:,7);
   vModelError = (mGamsStats(:,1) == 1) + (ismember(mGamsStats(:,2),[1 2 8])) ~= 2;
   vReturnFlag(vModelError) = -vReturnFlag(vModelError);

   %[vReturnFlag NaN(NumSolutions,1)  mLevelsOrd NaN(NumSolutions,1) mResultsInt]
   
   %[vObjsCurr, mResultsIntCurr, mResultsValCurr, uelsOut, vReturnFlagCurr, mGamsStatsCurr] = GamsPCShell(sGamsFile,vFixed,vFixedVal,vLToUse,uels,Tolerance,InputPassed,1);
    
    uMap = MapLabels(uelsOut,uels);
    
    %Map the decision variable results back out into the uels provided
    mResultsInt = mResultsInt(:,uMap);
    mResultsVal = mResultsVal(:,uMap);                       
end

