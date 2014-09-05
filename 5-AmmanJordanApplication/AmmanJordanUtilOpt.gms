$Title Utility Shortage Management Model
$inlinecom { }
$Ontext
   Solves two-stage probabilistic optimization program to identify optimal solution
   (lowest expected cost) mix of integrated
   short- and long-term water conservation and supply expansion management
   actions. Actions must meet a probablistic range of shortage events.

   Uses a deterministic model formulation (fixed data inputs and outputs, [Wilchfort and Lund, 1997])
   that minimizes the expected costs (minimize mean value)

   Numbers and actions for Amman, Jordan water utility circa 2004. Data collected during
   interviews, meetings, and follow-up visits in Amman, Jordan with 21 professionals
   working at MWI, WAJ, JVA, LEMA, USAID, and in private consultation from
   November 2005 through January 2006.

   Resolves model for uncertain demand data forecast for 2005, 2020, 2025, 2040 etc as specified by binary flags
   in the parameter YtoUse(y). Each year is indenpendent.

   MODEL INPUT DATA: AmmanJordanUtilData.gdx
   MODEL OUTPUTS: written to AmmanJordanUtilOpt.gdx

   FURTHER DESCRIPTION
   A full description of the data and model are available in:
   David E. Rosenberg and Jay Lund (2009). "Modeling Integrated Water Utility Decisions with Recourse and Uncertainty".
   Water Resources Management. 23 (1), pp. 85-115. doi: 10.1007/s11269-008-9266-4. http://www.springerlink.com/content/k7h71596u3065104/

   #####################
   Programmed by David E. Rosenberg

   Dept. of Civil & Env. Engineering and Utah Water Research Lab
   Utah State University
   david.rosenberg@usu.edu

   HISTORY
     - Programmed by David Rosenberg August 27, 2006 (U.C. Davis, Phd dissertation work)
     - Derives from UtilityTwoStateMYMV.gms
     - GDX data dump added February 2014
     - Read input from GDX file August 2014

   CITATION:
   David E. Rosenberg and Jay Lund (2009). "Modeling Integrated Water Utility Decisions with Recourse and Uncertainty".
   Water Resources Management. 23 (1), pp. 85-115. doi: 10.1007/s11269-008-9266-4. http://www.springerlink.com/content/k7h71596u3065104/

   LICENSING:
   This code is distributed AS-IS with no expressed or implied warranty regarding functionality. The entire code or parts
   may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires
   prior written permission from the author. Use requires downloading and installing the General Algebraic Modeling System (www.gams.com)
   as well as licenses for the BDMLP and DICOPT solvers to solve mixed-integer and non-linear mixed integer programs.

   BUG REPORTS and FEEDBACK:
   This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker on the
   GitHub repository where you downloaded this file. Note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.

$Offtext

OPTION limrow = 0
OPTION limcol = 0

***** Load the SETS at Execute time
$gdxin AmmanJordanUtilData
$LOAD
SETS  i Long term actions
      l(i) Long term conservation actions
      j Short term actions
      s Seasons
      e Hydrologic events
      m(j) conveyance actions for capacity expansion
      n(j) applied water for waste-water reuse
      sc scenario data
      y  demand year
      d  Data Uncertainty
      mvi Near optimal scenarios
      t time years
      cls Action classifications;
*      ModTyp Model Formulation Types
*                        /conm "constrain mean", conv "constrain variance", mobj "multi objective"/;
$LOAD i,l,j,s,e,m,n,sc,y,d,mvi,t,cls
$gdxin

Display i,l,j,s,e,m,n,sc,y,d,mvi,t;

Alias (e, e2);
Alias (s, s2);
Alias (i, i2);
Alias (j, j2);
*Alias (mvi, mvi2);

SCALARS
    tcap Current Treatment Plant Capacity (million m3 per year)
    disttype Distribution type (1=uniform  2=normal)
    r Discount rate
    SDWidth Width of data range (number of standard deviations)
    td Total demand (Mm3 per year)
    yd Year x demand (Mm3 per year)
    CThresh Cost threshold to end search (JD mill)
    CUST Number of customers
    mill One million /1000000/
    VLIMIT upper limit on variance
    MLIMIT upper limit on mean
    VHI   highest value of variance
    VLO   lowest value of variance
    MHI   highest value of mean
    MLO   lowest value of mean
    TRAD  tradeoff coefficient between mean and variance costs (fraction)
    GAMMA near optimal tolerance (fraction of optimal objective function value)
    FOPT  optimal value of objective function
    RunMode type of optimization problem to solve (1=optimum 2=bounds for unfixed variables 3=bounds for all variables 4=crawl all solutions) /1/
    InputPassed input passed to use as set values (1=levels [lFixedVal] 2=volumes [lFixedValVol]) /1/;

PARAMETERS
    cap(i,d) Capital Costs of long-term actions (JD Mill)
    life(i,d) Lifespan of long-term actions (years)
    Lmaxunits(i) Upper bound on long-term actions (units)
    ScWDem(i)    Upper limit scales with demand multiplier (0 or 1)
    Gadj(i,j) Interaction for long-term conservation action needs adjustment so short-term limit is greater than zero
    al(j)        Apparent loss (fraction of effectivess that is an apparent loss)
    ScWDemST(j)  Scales with demand for Short-term actions (fraction of demand)
    Lmaxp(i,d) Effectiveness of long-term actions (MCM per unit per year)
    C2(j,s,d) UCost of short-term actions (JD Mill per event)
    sf(i,s,d) Seasonal distribution of long-term action effectiveness (fraction)
    dem(e,s,d) Shortage in event e of season s (fraction of total demand)
    p(e,s,d) Probability of shortage in event e of season s (fraction)
    g(i,j,d) Interaction between long and short-term actions (fraction)
    Smaxp(j,s,d) Upper limit of short-term actions (MCM per event)
    rf(s,d)  return flow (fraction of applied water use)
    tl(s,d)  treatment loss (fraction of return flow)
    sd(s,d)  Seasonal distribution of demand (fraction of normal demand)
    demyr(y,e,d) Shortage in year y (fraction of total demand)
    pyr(y,e,d)     Probability of shortage event e in year y (fraction)
    basedem(y,d) Base demand in year y (MCM per year)
    OpCost(i,d) Operating cost for long-term actions that do not have short-term dummy (JD per m3)
    CstScWDem(i) Cost scales with demand (1 or 0)
    YtoUse(y)    Run model for demand year (1 = yes and 0 = no)
*    ModToUse(ModTyp)  Run tradeoff curve solution method (1 = yes and 0 = no)
    LFixed(i)    Long-term action is fixed (1 = fixed and 0=free)
    LTOUSE(i)    Long-term action to use in current near-optimal iteration (1=yes 0=no)
    LFixedVal(i)   Value to fix long-term action at (integer levels)
    LFixedValVol(i) Volume to fix long-term action at (MCM per year)
    LngCls(i,cls) Classification of long-term actions;


** To deal with variable number of data scenarios for robust model
Set bs Bases /b1 * b3/;

PARAMETER
   BaseVals(bs) /b1 1, b2 2, b3 5/
   ScenToUse(sc) Indicator of which scenarios to use;

SCALAR BaseIn Base index to start at
       BasePo Base power to start at;

** To deal with annualized discounting for long-term capital costs
PARAMETER AMORT(t) Amortization of capital cost;

Set CostTyp Types of cost /Long, Short, Total, Variance, MLimit, VLimit, SolStat, ModStat, NumScen, ItsUsd, NumEqu, NumVar, SolTime/
    Const   Constraint Type /Deliv2, Convey5, WWT6, prob, volum, variance/
    ModDat  Data for model /cost, level, volum, marg, maxvolum/
    Dum    Has one element /Val/;

********** Deterministic Model Formulation **********************
*     Deterministic versions of parameters, use average of Min and Max
      PARAMETERS
      cap_DET(i)       Capital cost of long-term actions (JD million)
      life_DET(i)      Lifespan of long-term actions (years)
      AMRTFC_DET(i)    Amortized fraction for long-term action (fraction)
      c1_DET(i)        Annualized capital cost for long-term actions (JD million per year)
      c2_DET(j,e,s)    Unit cost for short-term action (JD per m3 per event)
      p_DET(e,s)       Probability of shortage event e in season s (fraction)
      dem_DET(e,s)     Shortage in event e of season s (fraction of td)
      sf_DET(i,s)      Seasonal distribution of long-term action effectiveness (fraction)
      lmaxp_DET(i)     Effectivenes of long-term action (Mm3 per unit per year)
      smaxp_DET(j,e,s) Upper limit of short-term action (Mm3 per event)
      g_DET(i,j,s,e)       Supply enhancement interaction factor (1 or 0)
      rf_DET(s)        Return flow factor (fraction of applied water)
      tl_DET(s)        Treatment loss factor (fraction of return flow)
      sd_DET(s)        Seasonal distribution of demand (fraction of normal demand)
      ScDem_DET(e,s)   Demand scale factor (fraction above normal demand);

      Parameter SDef(s,e) Seasonal deficit
                SUse(s,e) Sesason to use;
      SCALAR    MinVal    Minimum value;


*     Error checking on upper limits for constraint 4
      Parameter E4_DET(j,e,s) Upper limit for deterministic equation 2
                E4_DS(j,e,s,sc) Upper limit for scenario equation 2;

*     Error checking on upper limits for constraint 2
      Parameter E2_DET(e,s) Upper limit for deterministic equation 2;

      VARIABLES
      DET_S(J,S,E)  Short term action for season S event E (million m3 per seas.)
      DET_L(I)      Long term action (units)
      DET_SLACK(S,E) Slack variable on meet shortage constraint (million m3 per seas.)
      DET_TCOST     Total expected annual cost for actions (JD Millions)
      DET_VCOST     Variance in annual cost for actions (JD millions)
      DET_MOCOST    Multi-objective combo of mean and variance costs (crazy units)
      L_VAL          Long term action value to use in current iteration of near optimal algorithm;

      POSITIVE VARIABLE   DET_S;
      INTEGER VARIABLE    DET_L;
*     Susinct way of representing non-negativity constraints [eq 7]

      EQUATIONS
      DEQ1        Objective function cost (JD Mill) [equation 1]
      DEQ2(S,E)   Action effectiveness must meet or exceed seasonal and event shortages (Mm3) [eq 2]
      DEQ3(I)     Long term actions less than limits (integer) [eq 3]
      DEQ4(J,S,E) Upper limits for short-term actions with supply enhancement or demand hardening interactions (Mm3 per season)[eq 4]
      DEQ5(S,E)   Distribution system capacity (Mm3 per season) [eq 5]
      DEQ6(S,E)   Mass balance on treated wastewater treated wastewater (Mm3 per season)[eq 6]
      DEQ7        Calculation of variance in cost (JD Mill) [eq 7]
      DEQ8        Constraint on variance in cost [eq 8]
      DEQ9        Constraint on mean cost [eq 9]
      DEQ10       Calculation of tradeoff in mean-variance costs [eq 10]
      DEQ11(S,E)  Action effectiveness plus slack variable must equal seasonal and event shortages (Mm3) [eq 11]
      DEQ12(S,E)  Either Slack variable or short term actions must be zero [eq 12]
      DEQ13       Near optimal tolerance (objective function units) [Eq 13]
      DEQ14       Long-term action value for current near optimal tolerance iteration (units of the long term action) [Eq 14];


      DEQ1..        SUM(I, c1_DET(I)*DET_L(I)) +  SUM((S,E,J), p_DET(E,S)*c2_DET(J,E,S)*DET_S(J,S,E)) =E= DET_TCOST;
      DEQ2(S,E)..   SUM(l, sf_DET(l,S)*(1+ScWDem(l)*ScDem_DET(e,s))*Lmaxp_DET(l)*DET_L(l)) + SUM(J, (1-al(j))*DET_S(J,S,E)) =G=
                      sd_DET(S)*td*dem_DET(E,S);
      DEQ3(I)..      DET_L(I) =L= Lmaxunits(I);
      DEQ4(J,S,E).. DET_S(J,S,E) =L= Smaxp_DET(J,E,S)*(1+ScWDemST(j)*ScDem_DET(e,s)) + SUM(I, sf_DET(i,s)*g_DET(i,j,s,e)*(1+ScWDem(i)*ScDem_DET(e,s))*Lmaxp_DET(i)*DET_L(i));
      DEQ5(S,E)..   SUM(m, DET_S(m,S,E)) =L= sf_DET("ExpandCap",S)*(TCAP + Lmaxp_DET("ExpandCap")*DET_L("ExpandCap"));
      DEQ6(S,E)..   DET_S("ReuseWW",S,E) =L= rf_DET(s)*tl_DET(s)*SUM(n, DET_S(n,S,E));

      DEQ7..        SUM((S,E), p_DET(E,S)*sum(J, abs(c2_DET(J,E,S)*DET_S(J,S,E)  -  SUM((S2,E2,J2), p_DET(E2,S2)*c2_DET(J2,E2,S2)*DET_S(J2,S2,E2)))**2)) =E=DET_VCOST;
      DEQ8..        SUM((S,E), p_DET(E,S)*sum(J, abs(c2_DET(J,E,S)*DET_S(J,S,E)  -  SUM((S2,E2,J2), p_DET(E2,S2)*c2_DET(J2,E2,S2)*DET_S(J2,S2,E2)))**2)) =L= VLIMIT;
      DEQ9..        SUM(I, c1_DET(I)*DET_L(I)) +  SUM((S,E,J), p_DET(E,S)*c2_DET(J,E,S)*DET_S(J,S,E)) =L= MLIMIT;
      DEQ10..       (1 - TRAD)*(SUM(I, c1_DET(I)*DET_L(I)) +  SUM((S,E,J), p_DET(E,S)*c2_DET(J,E,S)*DET_S(J,S,E))) +
                         (TRAD)*SUM((S,E), p_DET(E,S)*sum(J, abs(c2_DET(J,E,S)*DET_S(J,S,E)  -  SUM((S2,E2,J2), p_DET(E2,S2)*c2_DET(J2,E2,S2)*DET_S(J2,S2,E2)))**2))
                                  =E= DET_MOCOST;
      DEQ11(S,E)..   SUM(l, sf_DET(l,S)*(1+ScWDem(l)*ScDem_DET(e,s))*Lmaxp_DET(l)*DET_L(l)) + SUM(J, (1-al(j))*DET_S(J,S,E)) + DET_SLACK(S,E) =E=
                      sd_DET(S)*td*dem_DET(E,S);
      DEQ12(S,E)..  SUM(J, (1-al(j))*DET_S(J,S,E))* DET_SLACK(S,E) =E= 0;
      DEQ13..       DET_TCOST =L= GAMMA*FOPT;
      DEQ14..       SUM(i, LTOUSE(i)*DET_L(i)) =E= L_VAL;

* unconstrained model
      MODEL TWOSTAGE_DETUN /DEQ1, DEQ2, DEQ3, DEQ4, DEQ5, DEQ6, DEQ7, DEQ11, DEQ12/;
* near-optimal
      MODEL TWOSTAGE_DETNE /DEQ1, DEQ2, DEQ3, DEQ4, DEQ5, DEQ6, DEQ7, DEQ11, DEQ12, DEQ13, DEQ14/;

* min mean, constrain variance model
      MODEL TWOSTAGE_DETCO /DEQ1, DEQ2, DEQ3, DEQ4, DEQ5, DEQ6, DEQ7, DEQ8, DEQ11, DEQ12/;
* min variance, constrain mean model
      MODEL TWOSTAGE_DETMCO /DEQ1, DEQ2, DEQ3, DEQ4, DEQ5, DEQ6, DEQ7, DEQ9, DEQ11, DEQ12/;
* Multi-objective tradeoff model
      MODEL TWOSTAGE_DETMOBJ /DEQ1, DEQ2, DEQ3, DEQ4, DEQ5, DEQ6, DEQ7, DEQ10, DEQ11, DEQ12/;


*Solve Options
OPTION MIP = BDMLP;
OPTION MINLP = DICOPT;
*OPTION MINLP = BARON;
TWOSTAGE_DETUN.ITERLIM = 25000;
TWOSTAGE_DETNE.ITERLIM = 25000;
TWOSTAGE_DETCO.ITERLIM = 25000;
TWOSTAGE_DETMCO.ITERLIM = 25000;
TWOSTAGE_DETMOBJ.ITERLIM = 25000;

*TWOSTAGE_DS.ITERLIM = 50000;

******** DEFINITIONS OF SETS & PARAMETERS TO FORMAT OUTPUT ***********************
SET sol Solutions generated /sol1*sol10000/
    levs Level data /Level, LBnd, UBnd, CurrV, SolStat, ModStat, SolGen/;


*    use Dum as last parameter dimension to force into row normalized form
*    when output to Excel for use by pivot tables


PARAMETERS
    LongActs(y,sol,mvi,i,ModDat,Dum) Values for long-term actions
    MargLT(y,sol,mvi,i,Dum)  Reduced costs (marginals) for long-term actions
    Lives(i)   Lifespan of long-term actions
    ConstShad(y,Const,S,E,mvi,d,Dum) Shadow values on constraints
    ShortActs(y,j,s,e,mvi,d,dum) Data values for short-term actions in each data scenario
    BaseDemOut(y,d)      Base demand to output to excel
    Costs(y,sol,mvi,CostTyp,Dum) Cost parts for models
    SCCosts(y,d,mvi,CostTyp,dum) Value of cost component
    VLimValues(y,d,mvi) Upper limit on variance
    MLimValues(y,d,mvi) Upper limit on mean
    TradVals(y,mvi)     Tradeoff coefficient value
    LongActsOpt(y,i,ModDat,Dum) Optimal Values for long-term actions
    CostsOpt(y,CostTyp,Dum) Cost parts for Optimal solution
    LevelData(y,sol,mvi,Levs,Dum) Level data for iterations
    YtoUse(y) Data on enumeration levels
    STBound(j,s,e) Short term bound to use when a corresponding long-term action is fixed (MCM)
;

* LOAD in the all the prior results from the GDX file
Execute_load 'AmmanJordanUtilData.gdx',
    tcap
    disttype
    r
    SDWidth
    td
    yd
    CThresh
    CUST
    mill
    VLIMIT
    MLIMIT
    VHI
    VLO
    MHI
    MLO
    TRAD
    GAMMA
    FOPT
    cap
    life
    Lmaxunits
    ScWDem
    Gadj
    al
    ScWDemST
    Lmaxp
    C2
    sf
    dem
    p
    g
    Smaxp
    rf
    tl
    sd
    demyr
    pyr
    basedem
    OpCost
    CstScWDem
    YtoUse
    AMORT
    LngCls
    cap_DET
      life_DET
      AMRTFC_DET
      c1_DET
      c2_DET
      p_DET
      dem_DET
      sf_DET
      lmaxp_DET
      smaxp_DET
      g_DET
      rf_DET
      tl_DET
      sd_DET
      ScDem_DET
      SDef
      SUse
      MinVal
      E4_DET
      E4_DS;

*Loop over demand years
Loop (y$YtoUse(y),
*** Prepare input values

    p(e,s,d) = pyr(y,e,d);
    dem(e,s,d) = demyr(y,e,d);
    yd = (basedem(y,"min") + basedem(y,"max"))/2;

*   Deterministic input values
    dem_DET(e,s) = (dem(e,s,"max")+dem(e,s,"min"))/2;
    p_DET(e,s) = (p(e,s,"max")+p(e,s,"min"))/2;
    g_DET(i,j,s,e) = (g(i,j,"max")+g(i,j,"min"))/2;
    ScDem_DET(e,s) = yd / td - 1;
    ScDem_Det(e,s)$(yd le 0) = 0;
*   Adjust probability of last event or season so fractions sum to 1
    p_DET(e,s)$(ord(e) eq card(e)) = 1 - SUM(e2$(ord(e2) ne card(e2)),p_DET(e2,s));
*   Adjust interaction values of g to avoid infeasibilities (negatives right-hand sides for equation 4)
    Loop((i,j)$Gadj(i,j),
         SDef(s,e) = (Smaxp_DET(j,"e1",s)*(1+ScWDemST(j)*ScDem_DET(e,s)) + sum(i2$(g_DET(i2,j,s,e) lt 0),g_DET(i2,j,s,e)*(1+ScWDem(i2)*ScDem_DET(e,s))*Lmaxp_DET(i2)*sf_DET(i2,s)))/(Lmaxp_DET(i)*sf_DET(i,s)*(1+ScWDem(i)*ScDem_DET(e,s)));
         SUse(s,e) = 0;
         MinVal = 0;
         Display SDef;
*        Find minimum Sdef
         Loop ((s,e),
             If (SDef(s,e) lt MinVal,
                 MinVal = SDef(s,e);
                 );
             );
         If (MinVal lt 0,
            g_DET(i,j,s,e) = g_DET(i,j,s,e) - MinVal;
            );
         );

    Display g_det;
*   Set limits on short term variables

    DET_SLACK.LO(S,E) = -(SUM(l, sf_DET(l,S)*(1+ScWDem(l)*ScDem_DET(e,s))*Lmaxp_DET(l)*Lmaxunits(l)) -  sd_DET(S)*td*dem_DET(E,S));
    DET_SLACK.LO(S,E)$(DET_SLACK.LO(S,E) gt 0) = 0;
    DET_SLACK.UP(S,E) = 0;

*   Det_L.UP("ConsProg") = 0;
*   Det_L.LO("WWReuse") = 1;
*   Det_L.LO("DesalSW") = 1;
*    DET_L.UP("RedPhyLeak") = 1;
*    DET_L.LO("RedPhyLeak") = 1;

*     TradVals(y,mvi) =  1/(card(mvi) - 1)*(ord(mvi)- 1);
*Switch the second and last values to reflect ordering in Constrained Method mvi solution
*     TradVals(y,mvi)$(ord(mvi) gt 2) = TradVals(y,mvi-1);
*     TradVals(y,mvi)$(ord(mvi) eq 2) = 1;

*Display TradVals;

***      First solve for optimum
         DET_L.UP(I) = Lmaxunits(I);
         DET_L.L(I) = Lmaxunits(I);
         DET_SLACK.L(S,E) = 0;

         SOLVE TWOSTAGE_DETUN USING MINLP MINIMIZING DET_TCOST;

         FOPT = DET_TCOST.L;

***      Record optimal solution results in first i2 and d bins
         LongActsOpt(y,i,"level","Val") = DET_L.L(i);
         LongActsOpt(y,i,"cost","Val") = DET_L.L(i)*c1_DET(i);
         LongActsOpt(y,i,"marg","Val") = DET_L.M(i);
         LongActsOpt(y,i,"volum","Val") = (1+ScWDem(i)*Sum((e,s),ScDem_DET(e,s))/(card(e)*card(s)))*Lmaxp_DET(i)*DET_L.L(i);
         LongActsOpt(y,i,"maxvolum","Val") = (1+ScWDem(i)*Sum((e,s),ScDem_DET(e,s))/(card(e)*card(s)))*Lmaxp_DET(i)*Lmaxunits(i);
         CostsOpt(y,"Total","Val") = DET_TCOST.L;
         CostsOpt(y,"Long","Val") = SUM(I, C1_DET(I)*DET_L.L(I));
         CostsOpt(y,"Variance","Val") = DET_VCOST.L;
         CostsOpt(y,"SolStat","Val") = TWOSTAGE_DETUN.SOLVESTAT;
         CostsOpt(y,"ModStat","Val") = TWOSTAGE_DETUN.MODELSTAT;
         CostsOpt(y,"ItsUsd","Val") = TWOSTAGE_DETUN.ITERUSD;
         CostsOpt(y,"NumEqu","Val") = TWOSTAGE_DETUN.NUMEQU;
         CostsOpt(y,"NumVar","Val") = TWOSTAGE_DETUN.NUMVAR;
         CostsOpt(y,"SolTime","Val") = TWOSTAGE_DETUN.RESUSD;

         CostsOpt(y,"Short",dum) = CostsOpt(y,"Total",dum) - CostsOpt(y,"Long",dum);
         CostsOpt(y,"NumScen","Val") = 1;

         LongActsOpt(y,i,ModDat,dum) = LongActsOpt(y,i,ModDat,dum) + EPS;


    );
*** end of demand year loop

*Output input scalars*
Set scals Set whose elements are model scalars /tcap, disttype,r,sdwidth,td,numyrs,numscens, numltacts, numcsts, nummvpts, gamma, fopt/;

Parameter ScalVals(scals);

ScalVals("tcap") = TCAP;
ScalVals("disttype") = DISTTYPE;
ScalVals("r") = R;
ScalVals("sdwidth") = SDWIDTH;
ScalVals("td") = TD;
ScalVals("numyrs") = Sum(y,YtoUse(y));
ScalVals("numscens") = Card(sc);
*ScalVals("numscits") = 3 + Sum(ro$(NumScens(ro) gt 1), 1);
ScalVals("numltacts") = Card(i);
ScalVals("numcsts") = Card(CostTyp);
*ScalVals("CumScens") = Sum(ro$(NumScens(ro) gt 1), NumScens(ro));
ScalVals("nummvpts") = Card(mvi);
ScalVals("gamma") = GAMMA;
ScalVals("fopt") = FOPT;

*Dump results to GDX file
Execute_Unload 'AmmanJordanUtilOpt'
