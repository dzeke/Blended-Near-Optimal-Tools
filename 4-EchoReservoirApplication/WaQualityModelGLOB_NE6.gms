$ontext
Created by Omar Alminagorta
and supervised by Dr. David Rosenberg
Utah State University
May-2011

****** The Original Model ***********

This is a model to select the Best Management Practice(BMPs) based on
specific targets and considering economic and efficiency factors.
This model is applied at the Echo Reservoir and this is defined as follow:

- Objective function: Minimize Cost
- Decision Variable: Area to implement the BMP
- Constraints: Area or Lenght Available

A full description of the model is available at:

Omar Alminagorta, Bereket Tesfatsion, David E. Rosenberg, Bethany Neilson (2013). "Simple Optimization Method to Determine Best Management
Practices to Reduce Phosphorus Loading in Echo Reservoir, Utah". ASCE-Journal of Water Resources Planning and Management, 139(1), 122-125, doi: 10.1061/(ASCE)WR.1943-5452.0000224.
http://ascelibrary.org/doi/abs/10.1061/%28ASCE%29WR.1943-5452.0000224

***** Modifications to undertake Near-Optimal Analysis *****
David E. Rosenberg

1) Solve for optimal solution by several model formulation methods. Check that b, c, and d give the same solution

   a. Model WaterQualitySub : Original version with phosphorus removal targets specific to each sub-watershed [Scenario #1 in Alminagorta et al (2013)]
   b. Model WaterQuality : Original version with global phosphorus reduction target across sub-watersheds [Scenario #2 in Alminagorta et al (2013)]
   c. Model WaterQualityEx : Variant of WaterQuality that only includes decision variables allowed by the CAFF parameter (Used here) [This reduces the problem size
         from (10 BMPS)*(3 sources)*(3 sub-watersheds) = 90 down to sum(CAFF(i,s))*3 = 39 variables].
   d. Model WaterQualityFull : Variant of WaterQualityEx that uses full matrix of constraints (Ax <= b). This is a check as the full matrix of constraint values
        are later passed to Matlab and used in various near-optimal analyses
   e. Model WaterQualityRem : Variant of WaterQualityEx that includes a variable and equation for total phosphorus removed.
        this version is used to identify Pareto solutions between the Cost and Removal objectives by the constraint method

2) Solve for alternatives that are maximally-different in decision space (Modeling to Generate Alternatives (MGA) by the Hop-Skip-Jump approach).
         Brill, E. D., Jr., Chang, S.-Y., and Hopkins, L. D. (1982).
         "Modeling to Generate Alternatives: The HSJ Approach and an Illustration Using a Problem in Land Use Planning." Management Science, 28(3), 221-235.

    f. Model HSJEx : Adds Hop-Skip-Jump objective and near-optimal tolerance constraint to the prior WaterQualityEx model.
    g. Model MaxDiff: Adds maximally different objective and near-optimal tolerance constraint to the prior WaterQualityEx model.

3) Dump all the results to a .gdx file

4) These substitutions were made to the original problem to make computations easier
  a - Substitute ef(i)*Areaused(i,w,s) for P(i,w,s) to eliminate the efficiency constraint. Removes the variable P(i,w,s)
        from the model formulation.
  b - Relax REQUIRE constraint to allow trading among sub-watersheds.
  c - Explicitly add non-negativity constraints to allow movement along these edges
  d - Explictly excludes decisions that are not possible as specifiedi in the matrix Capp(i,s)

INPUTS
  All provided in this model file (no external)

OUTPUTS
  Results dumped to a gdx file (includes inputs as well)

* Version history

%% #####################
%   Programmed by David E. Rosenberg, Omar and Alminagorta
%
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
%
%   History
%    - 2011 - Omar Alminagorts developed the original model and scenarios a and b described above
     - January 18, 2013 -- David Roenberg added model versions c and d to make computations easier for subsequet near-optimal analysis
     - January 2014 -- David Rosenberg added model version f to compare results to the Modeling to Generate Altneratives method
     - July 2014 -- David Rosenberg added model version e to generate pareto optimal solutions by the constraint method to a multiobjective problem that
                 increases the total phosphorus removed.
%
%   Citation:
%   If using model versions a or b, cite the original work:
         Omar Alminagorta, Bereket Tesfatsion, David E. Rosenberg, Bethany Neilson (2013). "Simple Optimization Method to Determine Best Management
         Practices to Reduce Phosphorus Loading in Echo Reservoir, Utah". ASCE-Journal of Water Resources Planning and Management, 139(1), 122-125, doi: 10.1061/(ASCE)WR.1943-5452.0000224.
         http://ascelibrary.org/doi/abs/10.1061/%28ASCE%29WR.1943-5452.0000224

%   If using model versions c,d,e,f cite:
%        David E. Rosenberg (in review) "Near-optimal alternative generation,
%        visualization, and interaction for water resources decision making".
%        Water Resources Research. Submitted August 2014

%   Licensing:
%   This code and all the dependent functions are distributed AS-IS with
%   no expressed or implied warrenty regarding the claimed functionility. If you use this code, cite
%   according to the Citation information above.
%
%   Bug Reports and Feedback:
%   This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.
%
%% #######################

**** Description of the original model *****

g represent the available Resources
Agric Land(km2)               1,3,4,5
Grass filter strips (km2)       2
Grazing Land(km2)               6
Stream Bank(Km)                 7
Stream Fences (Km)              8
Manure System(System)           9


s1  Land Applied Manure
s2  Private Land Grazing
s3  Diffuse Run-off

w1  Chalk Creek
w2  Weber River below Wanship
w3  Weber River above Wanship

BMP
1        Land retirement
2        Grazing land protection
3        Stream fencing
4        Stream bank stabilization
5        Cover crops
6        Grass filter strips
7        Animal waste facility
8        Conservation tillage
9        Agricultural nutrient management
10        Sprinkler irrigation

More details see the manuscript approved in Jan2012: " Simple Optimization Method to Determine Best Management Practices to Reduce Phosphorus Loading in Echo Reservoir, Utah"
Alminagorta, O., Tesfatsion, B., Rosenberg, D., and Neilson, B.
http://ascelibrary.org/doi/abs/10.1061/%28ASCE%29WR.1943-5452.0000224
$offtext


set  i /BMP1  Retire land
        BMP2  Protect grazing land
        BMP3  Fence streams
        BMP4  Stabilize stream banks
        BMP5  Cover crops
        BMP6  Grass filter strips
        BMP7  Animal waste facility
        BMP8  Conservation tillage
        BMP9  Manage ag. nutrients
        BMP10 Sprinkler irrigation/

*   i BMPs / BMP1*BMP10/
    s Source of N / s1 "Manure", s2 "Grazing", s3  "Runoff"/
*    w Subwatershed  /w1*w3/
    w Subwatershed  /w1 "Chalk Creek", w2 "Weber below Wanship", w3 "Weber above Wanship"/
    g row to generate the inplementation matrix /1*9/
    ctype constraint type /request, existload/
;

ALIAS (i, i2);
ALIAS (s, ss2);
ALIAS (w, w2);

SCALAR GAMMA Near Optimal Tolerance (fraction of optimal objective function value) /1.10/
       NEb Near optimal limit (units of objective function value);

Table target(w,s)

              S1            S2                s3
w1          768.80        353.50            915.46
w2          753.60        155.00            549.30
w3         2848.00        371.80           1351.70

*Existing Total Watershed Load (kg/year)

Table ExisLoad(w,s)

          S1          S2          s3
w1        961        3535        6539
w2        942        1550        5493
w3        3560       3718        13517

;

Table area_availa (g,w)
*In the row 1: means that in subwatershed1 (W1), There are only 15.6Km2 of Available Agricultural Area

          W1          W2          W3
1        15.6        24.5        69.8
2        3.5         1.8         3.5
3        15.6        24.5        69.8
4        15.6        24.5        69.8
5        15.6        24.5        69.8
6       311.8        95.3       229.4
7        70.0        36.0        70.0
8        70.0        36.0        70.0
9         0.0         0.0         0.0

;

Parameter ef(i) reduction of P in Kg per year per unit of action implement (kg per area or kg per stream bank length)
/
BMP1        184.94
BMP2        21.74
BMP3        23.81
BMP4        208.34
BMP5        43.71
BMP6        143.47
BMP7        27106.66
BMP8        20.18
BMP9        43.71
BMP10       2914.21
/;

Parameter cost(i) cost of implementing action i ($ per kg phosphorus removed)

/
BMP1        7367.85
BMP2        238.10
BMP3        1373.48
BMP4        15.43
BMP5        674.61
BMP6        412.26
BMP7        729.73
BMP8        337.31
BMP9        167.55
BMP10       127.87
/
;

Variables

Areaused(i,w,s) area to be used by BMP
REMOVE total phosphorus removed (kg)
*P(i,w,s) phosporus removed
obj
;
*Positive Variables
*Areaused, P;
*Areaused;

VARIABLE PhosRemBounds Bounds on phosphorus removal used for pareto optimal analysis;

table E(g,i)  Exclusion matrix for row

        BMP1    BMP2    BMP3    BMP4     BMP5   BMP6   BMP7    BMP8   BMP9   BMP10
1        1                                1
2                                                1
3        1                                                      1
4        1                                                              1
5        1                                                                     1
6                1
7                        1
8                                1
9                                                        1
;

table Capp(i,s)  Binary parameter 1 if BMP can be applied to source s

            s1      s2       s3
BMP1        0        0        1
BMP2        0        1        0
BMP3        0        1        0
BMP4        0        0        1
BMP5        0        0        1
BMP6        1        1        1
BMP7        0        0        0
BMP8        1        0        1
BMP9        1        0        1
BMP10       0        0        1


;


Equations
*Effic(i,w,s) //Substitute definition of P to eliminate
**Request(w,s) Meet removal target
Request(s) Meet source targets across sub-watersheds - global version
RequestSub(w,s) Meet source targets in individual sub-watersheds
Are(g,w) Available resources
ExistingLoad(w,s) Stay at or below existing load
NonNeg(w,s,i) Decision variables must be positive
Objective Cost minimizing objective function
**Excluded versions of the equations
RemovedEx Total phosphorus removed by all BMPs sources and sub-watersheds
RequestEx(s) Meet source targets across sub-watersheds - global version
RequestSubEx(w,s) Meet source targets in individual sub-watersheds
AreEx(g,w) Available resources
ExistingLoadEx(w,s) Stay at or below existing load
NonNegEx(w,s,i) Decision variables must be positive
ObjectiveEx Cost minimizing objective function
;

*Quantity of Phosphorus removed by BMP in Kg/yr
*Effic(i,w,s)..   P(i,w,s) =e=ef(i)*Areaused(i,w,s);

* Meet Removal Target


Request(s)..
        sum((i,w),ef(i)*Areaused(i,w,s)*Capp(i,s)) =g= sum(w,target(w,s));

RequestSub(w,s)..
     sum(i,ef(i)*Areaused(i,w,s)*Capp(i,s)) =g= target(w,s);

* Available Resources
Are(g,w)..
         sum((s,i),(Capp(i,s)*E(g,i)*Areaused(i,w,s)))=l= area_availa(g,w);


*constraint to get result lower than existing load  Capp(i,s) represented the source matrix
ExistingLoad(w,s)..

          sum((i),ef(i)*Areaused(i,w,s)*Capp(i,s)) =l= ExisLoad(w,s);

*Non-negativity constraints
NonNeg(w,s,i).. Areaused(i,w,s) =g= 0;

Objective..  obj=e= sum((i,w,s),ef(i)*Areaused(i,w,s)*cost(i));

*****
***** Excluded versions of the equations

RemovedEx..
        REMOVE =e= sum((i,w,s)$CAPP(i,s),ef(i)*Areaused(i,w,s));


RequestEx(s)..
        sum((i,w)$CAPP(i,s),ef(i)*Areaused(i,w,s)) =g= sum(w,target(w,s));

RequestSubEx(w,s)..
     sum(i$Capp(i,s),ef(i)*Areaused(i,w,s)) =g= target(w,s);

* Available Resources
AreEx(g,w)..
         sum((s,i)$Capp(i,s),(E(g,i)*Areaused(i,w,s)))=l= area_availa(g,w);


*constraint to get result lower than existing load  Capp(i,s) represented the source matrix
ExistingLoadEx(w,s)..

          sum((i)$Capp(i,s),ef(i)*Areaused(i,w,s)) =l= ExisLoad(w,s);

*Non-negativity constraints
NonNegEx(w,s,i)$Capp(i,s).. Areaused(i,w,s) =g= 0;

ObjectiveEx..  obj=e= sum((i,w,s)$Capp(i,s),ef(i)*Areaused(i,w,s)*cost(i));



option limrow = 140 ;
Model WaterQuality Original version with global removal target /Request,Are,ExistingLoad,NonNeg,Objective/ ;
Model WaterQualityEx Global version but excludes decisions not allowed by CAPP /RequestEx,AreEx,ExistingLoadEx,NonNegEx,ObjectiveEx/ ;
Model WaterQualityRem  Global exclude version with total removal /RemovedEx, RequestEx,AreEx,ExistingLoadEx,NonNegEx,ObjectiveEx/ ;
Model WaterQualitySub Original version with sub-watershed removal targets /RequestSub,Are,ExistingLoad,NonNeg,Objective/ ;

OPTIONS ITERLIM=10000;
OPTIONS RESLIM=1000000;

* This allows comments to be made in a line by using the characters: {}
$inlinecom{ }

* The maximum number of characters recognized in a single line is 400
$maxcol 400

*Allow empty sets if no nodes (intermediary points) defined

SETS
  j constraints /j1*j79/
  v vertices /v1*v250/
  rt results type /obj, modstat, solstat, vgen, jact, TotRem, Dist/
  ct constraint type groups /ct1*ct5/;

ALIAS(j, j2);
ALIAS(v, v2);
ALIAS(ct, ct2);

PARAMETERS
    VERTS(v,i,w,s) Coordinates of identified verticies
    BINDS(v,j) Binding constraints associated with vertex v (1=bind =not)
    BINDS1(v,s) Binding request associated with vertex v (1=bind 0=not)
    BINDS1B(v,w,s) Binding existing load associated with vertex v (1=bind 0=not)
    BINDS2(v,g,w) Binding are constraints associated with vertex v (1=bind 0=not)
    BINDS3(v) Binding near optimal constraints associated with vertex v (1=bind 0=not)
    BINDS4(v,w,s,i) Binding non-negativity constraint with vertex v (1=bind 0=not)
    CBINDS(j)  Binding constraints for current vertex minus current constraint (1=bind 0=not)
    CONTEST(j) Model status for constraint test results
    ObjVal(v,rt)  Model results at the vertex
    CTCOUNT(ct) Counting of constraints in each group (number)
    CTCUM(ct)   Cumulative count of constraints in current and prior groups (number)
    CUMCAPP(w,s,i) Cumulative count of BMPS that can be applied to source s (number) - used in indexing constraints
    CUMCAPPS(s) Number of BMPS that can be applied to source s (number)
    NonZeroDec(v) Number of non-zero decision variables at the vertex
    c(j)   Dummy cost coefficients associated with constraints (1 or 0)
    c1(s) Dummy cost coefficients associated with request constraints (1 or 0)
    c1B(w,s) Dummy cost coefficients associated with exisiting load constraints (1 or 0)
    c2(g,w) Dummy cost coefficients associated with are constraints (1 or 0)
    c3 Dummy cost coefficient associated with near optimal constraint (1 or 0)
    c4(w,s,i) Dummy cost coefficients associated with non-negativity constraint (1 or 0)
    AMat(j,w,s,i) Matrix of constraint coefficients (j constraints by (w x s x i) decisions))
    vB(j) right hand side constraints;

SCALARS
    cV Current vertex /1/
    lV Last vertex added /1/
    cJ Current constraint /1/
    NumSolvs Number of Solves /1/
    NumDims Number of dimensions (number of decision variables)
    NumPermutedDims Number of permuted decision variables (multiplying all indexes together)
    MaxVerts Maximum number of vertices

    MAXNUM very large number /2000000/;

***** Now to the Hop-Skip-Jump MGA approach to identifying maximally different near-optimal solutions
* Maximal different means minimize the sum of non-zero variables in the prior alternative

PARAMETERs NONZEROVAR(i,w,s) Decision variables that are non zero (1=non zero 0=zero)
           vNONZEROVAR(v,i,w,s) Decision variables that are non zero by vertex (1=non zero 0=zero)
           NumHSJSols Number of Hop-Skip-Jump alternatives generated (#)
           HSJStart Solution number to begin Hop-Skip-Jump iterations (#);
VARIABLES  HSJVAL Value of the hop-skip-jump objective function (decision variable units)
           S3 Slack and surplus variable for near optimal constraint (constraint units)

POSITIVE VARIABLES S3;

EQUATIONS
  HSJObj  Hop-Skip-Jump objective function
  NearOptimalCons Near optimal constraint;

HSJObj.. HSJVAL =E= sum((i,w,s)$(NONZEROVAR(i,w,s) and Capp(i,s)),AreaUsed(i,w,s));
NearOptimalCons.. Obj + S3 =E= NEb;

MODEL HSJEx Hop Skip Jump MGA /RequestEx,AreEx,ExistingLoadEx,NonNegEx,HSJObj,NearOptimalCons,ObjectiveEx/ ;
Model WaterQualityRemNE  Global exclude version with total removal and near-optimal tolerance constraint /RemovedEx, RequestEx,AreEx,ExistingLoadEx,NonNegEx,ObjectiveEx,NearOptimalCons/ ;


***** Now to the MaxDifferent approach to identifying maximally different near-optimal solutions
* Maximal different means maximize the distance between the current alternative and all other alternatives.
* Distance is measured as the sum of the absoluate values of the differences between the decision variable values of the new alternative and all prior alternatives
* The first formulation solves as a non-linear program using abs( )

SET        dir Direction /neg, pos/;
SCALAR     DistCriteria  Stopping distance for iterations (kg phosphorus changed) /100/
           CurrDist Distance objective function value of current iternation (kg phosphorus changed)
           UpperBound Upper bound on the absolute value of the difference (kg phosphorus changed);
PARAMETERS  VertsToUse(v) Verticies to include in current iteration (1=yes 0=no)
            AllDists(v,v2) Distance between vert v and v2 (kg phosphorus removed)
            DirVal(dir) Directional multiplier /neg -1, pos 1/
            ContMGA Binary values to indicate to continue MGA iterations (1=continue 0=stop)
           UpperBoundVar(i,w,s) Upper bound of the specified variable (kg phosphorus changed)
           VarToUse(i,w,s) Binary indicator of variable to use (1=yes 0=no)
;

VARIABLES  DistVal Value of the distance objective function (kg phosphorus changed)
           Dist(v) Summed distance among components between current alternative and alternative v (kg phosphorus removed)
           DistComp(v,i,w,s) Component of distance between current alternative and alternative v on decision variable i w s (kg phosphorus changed)
           LevVal Level of selected decision variable (kg phosphorus removed)
Positive Variables Dist;

EQUATIONS
  OverallDiff(v) Overall distance is L1 minimum distance between alternatives
  AbsDiff(v) Absolute value of difference
  DistBnds(v) Upper bound on absolute value of difference
  DistCompDeff(v) Definition of distance components between current alternative and alternative v
  CurrVar Level of decision variable
      ;

OverallDiff(v)$VertsToUse(v)..   DistVal =L= Dist(v);
AbsDiff(v)$VertsToUse(v)..      Dist(v) =L= sum((i,w,s)$Capp(i,s),ef(i)*abs(AreaUsed(i,w,s)-VERTS(v,i,w,s)));
DistBnds(v)$VertsToUse(v)..      Dist(v) =L= UpperBound;
CurrVar..                      LevVal =E= sum((i,w,s)$Capp(i,s),VarToUse(i,w,s)*ef(i)*AreaUsed(i,w,s));

MODEL VarBounds Find bounds for a specified variable /RequestEx,AreEx,ExistingLoadEx,NonNegEx,CurrVar,NearOptimalCons,ObjectiveEx/ ;
MODEL DistExAbs Distance MGA using absolute value /RequestEx,AreEx,ExistingLoadEx,NonNegEx,OverallDiff,AbsDiff,DistBnds,NearOptimalCons,ObjectiveEx/ ;

*** A second formulation solves as a mixed integer program with DistCompDir representing the distance component in the positive and negative directions
* and DistBin the associated binary variable for the component

PARAMETERS DirOffset(dir) Offset to calculate upper bound on DistCompDir /pos 0, neg 1/;
VARIABLES Y(v,i,w,s) signed distance between current alternative and alternative v on decision variable i w s (kg phosphorus change from v)
          DistCompDir(v,i,w,s,dir) distance component in the specified direction (kg phosphorus change)
          DistBin(v,i,w,s) binary variable representing DistCompDir is in the positive direction (1=yes 0=no)

Positive Variables DistCompDir;
Binary Variables DistBin;

EQUATIONS
    BinDiff(v) Difference using binary approach
    SignedDiff(v,i,w,s) Signed value of difference
    DistBndsBin(v,i,w,s,dir) Upper Bound on differences
    DistDeff(v,i,w,s) Definition of the distance between alternative components  ;

BinDiff(v)$VertsToUse(v)..      Dist(v) =L= sum((i,w,s,dir)$Capp(i,s),DistCompDir(v,i,w,s,dir));
SignedDiff(v,i,w,s)$(VertsToUse(v) and Capp(i,s))..  Y(v,i,w,s) =E= ef(i)*(AreaUsed(i,w,s) - VERTS(v,i,w,s));
DistBndsBin(v,i,w,s,dir)$(VertsToUse(v) and Capp(i,s))..  DistCompDir(v,i,w,s,dir) =L=  UpperBoundVar(i,w,s)*(DirVal(dir)*DistBin(v,i,w,s) + DirOffset(dir));
DistDeff(v,i,w,s)$(VertsToUse(v) and Capp(i,s))..  sum(dir,DirVal(dir)*DistCompDir(v,i,w,s,dir)) =E= Y(v,i,w,s) ;

MODEL DistExBin Distance MGA using binary variables /RequestEx,AreEx,ExistingLoadEx,NonNegEx,OverallDiff,BinDiff,DistBndsBin,DistDeff,SignedDiff,NearOptimalCons,ObjectiveEx/ ;



* Initialize everthing
VERTS(v,i,w,s) = 0;
BINDS(v,j) = 0;
BINDS1(v,s) = 0;
BINDS1B(v,w,s) = 0;
BINDS2(v,g,w) = 0;
BINDS3(v) = 0;
BINDS4(v,w,s,i) = 0;

**** Formulate the overarching matrix describing the system of linear equations
*** AMat*Areaused <= vB
*** where AMat is an j x w x s x i matrix, Areaused is a w x s x i vector, and vB is j x 1 vector

**First develop the constraint counts
* S Request
CTCOUNT("ct1") = card(s);
*Existing Load
CTCOUNT("ct2") = card(s)*card(w);
*Are available resources
CTCOUNT("ct3") = card(g)*card(s);
*Near optimal
CTCOUNT("ct4") = 1;
*Non-negativity
CTCOUNT("ct5") = card(w)*sum((i,s),Capp(i,s));

CTCUM(ct) = sum(ct2$(ord(ct2) le ord(ct)),CTCOUNT(ct2));

NumDims = sum((i,s),Capp(i,s));
CUMCAPPS(s) = SUM(i,CAPP(i,s));

CUMCAPP(w,s,i)$CAPP(i,s) = (ord(w)-1)*NumDims + sum(ss2$(ord(ss2) lt ord(s)),CUMCAPPS(ss2)) + SUM((i2)$((ord(i2) le ord(i))),CAPP(i2,s));

Display CTCOUNT, CTCUM, CAPP, NumDims, CUMCAPP;

NumDims = card(w)*NumDims;

* 1. Work on the S Request Constraints (minus sign reverses direction of constraint)
AMat(j,w,s,i)$((ord(j) le card(s)) and (ord(j) eq ord(s)) and CAPP(i,s)) = -ef(i);
vB(j)$(ord(j) le card(s)) = - sum((w,s)$(ord(s) eq ord(j)),target(w,s));

* 2. Work on the Existing Load constraints
AMat(j,w,s,i)$((ord(j) gt card(s)) and (ord(j) le CTCUM("ct2")) and (ord(j) eq card(s)+(ord(w)-1)*card(s)+ord(s)) and CAPP(i,s)) = ef(i);
vB(j)$((ord(j) gt card(s)) and (ord(j) le CTCUM("ct2"))) = sum((w,s)$(ord(j) eq card(s)+(ord(w)-1)*card(s)+ord(s)),ExisLoad(w,s));

* 3. Work on the Are Available Resource Constraints
AMat(j,w,s,i)$((ord(j) gt CTCUM("ct2")) and (ord(j) le CTCUM("ct3")) and Capp(i,s)) = sum((g)$(ord(j)-card(s)-card(w)*card(s) eq (ord(g)-1)*card(w)+ord(w)),E(g,i));
vB(j)$((ord(j) gt CTCUM("ct2")) and (ord(j) le CTCUM("ct3"))) = sum((g,w)$(ord(j)-card(s)-card(w)*card(s) eq (ord(g)-1)*card(w)+ord(w)),area_availa(g,w));

* 4. Work on the near optimal constraint
AMat(j,w,s,i)$((ord(j) eq CTCUM("ct4")) and Capp(i,s)) = ef(i)*cost(i);
*Need to change after solving
vB(j)$(ord(j) eq CTCUM("ct4")) = -10;

* 5. Work on the non-negativity constraint (this should just be the negative of the idendity matrix
AMat(j,w,s,i)$((ord(j) gt CTCUM("ct4")) and (ord(j)-CTCUM("ct4") eq CUMCAPP(w,s,i)))= -1;
vB(j)$(ord(j) gt CTCUM("ct4")) = 0;

Display AMat, vB;

EQUATIONS
    FullCons(j) Full constraint set;

FullCons(j)..
     sum((i,w,s)$CAPP(i,s),AMat(j,w,s,i)*Areaused(i,w,s)) =L= vB(j);

Model WaterQualityFull /FullCons,ObjectiveEx/;
MODEL HSJFull Hop SKip Jump MGA Full version /HSJObj,ObjectiveEx, FullCons/;

**** Now to actually running ****

**** STEP 1 ********
**** Solve for optimal using the exclusions model version
Solve WaterQualityEx USING LP Minimizing obj;


**** record stuff about the solution
VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
ObjVal(v,"modstat")$(ord(v) eq cV) = WaterQualityEx.Modelstat;
ObjVal(v,"solstat")$(ord(v) eq cV) = WaterQualityEx.Solvestat;
ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));
*Store the phoshorus removal amount as lower bound for later use in pareto-optimal analysis
PhosRemBounds.LO = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));

*Determine which decision variables are basic (non-zero)
NONZEROVAR(i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)) = 1;
vNONZEROVAR(v,i,w,s)$(ord(v) eq cV) = NONZEROVAR(i,w,s);

NumPermutedDims = card(i)*card(s)*card(w);
NonZeroDec(v)$(ord(v) eq cV) = sum((i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)),1);
MaxVerts = Card(v);

* Establish right hand side for near-optimal tolerance constraint
NEb = GAMMA*sum(v$(ord(v) eq cV),ObjVal(v,"obj"));
vB(j)$(ord(j) eq CTCUM("ct4")) = NEb;



DISPLAY NEb, vB;

* Solve the original version problem
cV = cV+1;
lV = lV+1;
Solve WaterQuality USING LP Minimizing obj;

**** record stuff about the solution
VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
ObjVal(v,"modstat")$(ord(v) eq cV) = WaterQuality.Modelstat;
ObjVal(v,"solstat")$(ord(v) eq cV) = WaterQuality.Solvestat;
ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));

*Determine which decision variables are basic (non-zero)
NONZEROVAR(i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)) = 1;
NonZeroDec(v)$(ord(v) eq cV) = sum((i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)),1);
vNONZEROVAR(v,i,w,s)$(ord(v) eq cV) = NONZEROVAR(i,w,s);

*DISPLAY Are.SLACK;
CBINDS(j)=0;
*Map the separate constraints into the single j index

CBINDS(j)$((ord(j) le card(s)) and (sum((s)$(ord(j) eq ord(s)),Request.SLACK(s)) eq 0)) = 1;
*DISPLAY CBINDS;
CBINDS(j)$((ord(j) gt card(s)) and (ord(j) le CTCUM("ct2")) and (sum((w,s)$(ord(j)-card(s) eq (ord(w)-1)*card(s)+ord(s)),ExistingLoad.SLACK(w,s)) eq 0)) = 1;
*DISPLAY CBINDS;
CBINDS(j)$((ord(j) gt CTCUM("ct2")) and (ord(j) le  CTCUM("ct3")) and (sum((g,w)$(ord(j)-(card(s) + card(s)*card(w)) eq (ord(g)-1)*card(w)+ord(w)),Are.SLACK(g,w)) eq 0)) = 1;
CBINDS(j)$((ord(j) gt CTCUM("ct4")) and (ord(j) le  CTCUM("ct5")) and (sum((w,s,i)$(ord(j)-CTCUM("ct4") eq CUMCAPP(w,s,i)),NonNeg.SLACK(w,s,i)) eq 0)) = 1;
DISPLAY CBINDS ;
*CBINDS(j)$((ord(j) eq card(j)) and (NearOptimalCons.SLACK eq 0)) = 1;

* Solve the matrix version problem to identify to verify the matrix calcs are correct
cV = cV+1;
lV = lV+1;
Solve WaterQualityFull USING LP Minimizing obj;

**** record stuff about the solution
VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
ObjVal(v,"modstat")$(ord(v) eq cV) = WaterQualityFull.Modelstat;
ObjVal(v,"solstat")$(ord(v) eq cV) = WaterQualityFull.Solvestat;
ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));

*Determine which decision variables are basic (non-zero)
NONZEROVAR(i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)) = 1;
NonZeroDec(v)$(ord(v) eq cV) = sum((i,w,s)$((Areaused.L(i,w,s) ne 0) and CAPP(i,s)),1);
vNONZEROVAR(v,i,w,s)$(ord(v) eq cV) = NONZEROVAR(i,w,s);

**** STEP 2 ****
**** Solve by an MGA method
PARAMETER MGAMethod MGA Method (1=Hop Skip Jump  2=Max Difference NLP Abs 3=Max Difference MIP) /1/
         ExecTime       Execution time
         ParetoSols Number of pareto-optimal alternatives /10/;

*Start the time counter
ExecTime = timeExec;

if (MGAMethod eq 1,
*   *** Hop-Skip-Jump iterate to generate alternatives until one of three stopping criteria is reached
*   1) # basic variables = total variables
*   2) # basic variables on current iteration = # of basic variables on prior iteration
*   3) Iterations equal max. number of iterations (MaxVerts)

*   We will use the index v representing vertices (indices) as our iteration
*   At each iteration, recalculate the set of basic variables as all basic (non zero) variables
*   in the current alternative and prior alternative

     HSJStart = cV;

     while (((cV eq HSJStart) or ((cV le MaxVerts) and (sum(v$(ord(v) eq cV),NonZeroDec(v)) lt NumDims) and (sum(v$(ord(v) eq cV),NonZeroDec(v)-NonZeroDec(v-1)) gt 0))),

*    Increment the solution counters
          cV = cV+1;
          lV = lV+1;

           Solve HSJEx USING LP Minimizing HSJVAL;
*           Solve HSJFull USING LP Minimizing HSJVAL;

*   ****  record stuff about the solution
           VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
           ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
           ObjVal(v,"modstat")$(ord(v) eq cV) = HSJEx.Modelstat;
           ObjVal(v,"solstat")$(ord(v) eq cV) = HSJEx.Solvestat;
           ObjVal(v,"vgen")$(ord(v) eq cV) = HSJVAL.L;
           ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));

*     **  Calculate the nonzero decision variables to use in the next iteration
**        Add to the non-zero set variables that were previously zero but have non-zero values in this iteration.
*         Decision variables that were previously non-zero stay in the non-zero set even if
*         their value in the current solution is zero
           NonZeroVar(i,w,s)$((NonZeroVar(i,w,s) eq 0) and (Areaused.L(i,w,s) ne 0)) = 1;
           NonZeroDec(v)$(ord(v) eq cV) = sum((i,w,s),NonZeroVar(i,w,s));
           vNONZEROVAR(v,i,w,s)$(ord(v) eq cV) = NONZEROVAR(i,w,s);
           );
**         End Hop-Skip-Jump iterations

     NumHSJSols = cV - HSJStart ;
     );

if (MGAMethod gt 1,
*** Use maximize distance objective to iterate to generate alternatives until a stopping criteria is reached
*   1) Distance falls below some distance threshhold
*   2) Iterations equal max. number of iterations

*   We will use the index v representing vertices (indices) as our iterations. cV represents the current iteration
*   and HSJStart-1 to cV-1 the set of prior alternatives.
*   At each iteration, recalculate the set of prior alternatives from which to calculate distances

*   First some preprocessing to find the upper bound (maximum extent) of each decision variable
    LOOP((i2,w2,ss2)$CAPP(i2,ss2),
        VarToUse(i,w,s) = 0;
        VarToUse(i2,w2,ss2) = 1;

        SOLVE VarBounds USING LP MAXIMIZING LevVal;
*       Log the bound for this variable
        UpperBoundVar(i2,w2,ss2) = LevVal.L;
        );

*    Record the overall upper bound
     UpperBound=sum((i,w,s),UpperBoundVar(i,w,s));

     DISPLAY UpperBound,UpperBoundVar;

     HSJStart = cV;
     VertsToUse(v) = 0;
     ContMGA = 1;


*    Second iteratively solve the appropriate MGA model
     if (MGAMethod eq 2,
           DistExAbs.SolPrint = 1;
*          Run the MGA iterations solving the NLP version of the absolute value of the distance problem
           while (ContMGA ne 0,

*              Define the set of prior alternatives
               VertsToUse(v)$(ord(v) eq cV) = 1;
*              Increment the solution counters
               cV = cV+1;
               lV = lV+1;

*              Initialize solution at geometric mean of prior alternatives
               AreaUsed.L(i,w,s)$Capp(i,s) = sum(v$VertsToUse(v),Verts(v,i,w,s))/sum(v,VertsToUse(v));

               Solve DistExAbs USING DNLP Maximizing DistVAL;

****           record stuff about the solution
               VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
               ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
               ObjVal(v,"modstat")$(ord(v) eq cV) = DistExAbs.Modelstat;
               ObjVal(v,"solstat")$(ord(v) eq cV) = DistExAbs.Solvestat;
               ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
               ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));
               ObjVal(v,"Dist")$(ord(v) eq cV) = DistVal.L;
               AllDists(v,v2)$((ord(v) eq cV) and VertsToUse(v2)) = Dist.L(v2);
               ContMGA = (cV lt MaxVerts-ParetoSols) and (DistVal.L ge DistCriteria) and ((DistExAbs.Modelstat eq 1) or (DistExAbs.Modelstat eq 2) or (DistExAbs.Modelstat eq 7)) and ((DistExAbs.Solvestat eq 1) or (DistExAbs.Solvestat eq 2) or (DistExAbs.Solvestat eq 4));
               );
          NumHSJSols = cV - HSJStart;
          );

     if (MGAMethod eq 3,
*        Run the MGA iterations solving the MIP version of the absolute value of the distance problem
         DistExBin.SolPrint = 1;
         while (ContMGA ne 0,

*              Define the set of prior alternatives
               VertsToUse(v)$(ord(v) eq cV) = 1;
*              Increment the solution counters
               cV = cV+1;
               lV = lV+1;

               Solve DistExBin USING MIP Maximizing DistVAL;

****           record stuff about the solution
               VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
               ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
               ObjVal(v,"modstat")$(ord(v) eq cV) = DistExBin.Modelstat;
               ObjVal(v,"solstat")$(ord(v) eq cV) = DistExBin.Solvestat;
               ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
               ObjVal(v,"TotRem")$(ord(v) eq cV) = sum((i,w,s)$CAPP(i,s),ef(i)*Areaused.L(i,w,s));
               ObjVal(v,"Dist")$(ord(v) eq cV) = DistVal.L;
               AllDists(v,v2)$((ord(v) eq cV) and VertsToUse(v2)) = Dist.L(v2);
               ContMGA = (cV lt MaxVerts-ParetoSols) and (DistVal.L ge DistCriteria) and ((DistExBin.Modelstat eq 1) or (DistExBin.Modelstat eq 2) or (DistExBin.Modelstat eq 8)) and ((DistExBin.Solvestat eq 1) or (DistExBin.Solvestat eq 2) or (DistExBin.Solvestat eq 4));
               );

         NumHSJSols = cV - HSJStart-1;
         );
     );
** End Maximize Difference iterations
ExecTime = timeExec - ExecTime;


** STEP 3 ****
** Multi-objective analysis

*** First maximize the second objective to max Phosphorus remove
* Increment the solution counters
cV = cV+1;
lV = lV+1;

*Solve WaterQualityRem USING LP Maximizing REMOVE;
Solve WaterQualityRemNE USING LP Maximizing REMOVE;

**** record stuff about the solution
VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
ObjVal(v,"modstat")$(ord(v) eq cV) = WaterQualityRemNE.Modelstat;
ObjVal(v,"solstat")$(ord(v) eq cV) = WaterQualityRemNE.Solvestat;
ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
ObjVal(v,"TotRem")$(ord(v) eq cV) = REMOVE.L;
*Store the phoshorus removal amount as upper bound for later use in pareto-optimal analysis
PhosRemBounds.UP = REMOVE.L;

*** Use the constraint-based method to carve out the specified number of pareto-optimal alternatives between the PhosRemBounds
lV = 0;
while ((lV lt ParetoSols) and (cV le MaxVerts),
*   lV counts the iterations, cV counts the vertices and includes prior solves
*   increment the counters
    lV = lV+1;
    cV = cV+1;

*   Fixe remove to the current level
    REMOVE.UP = PhosRemBounds.UP - lV*(PhosRemBounds.UP - PhosRemBounds.LO)/ParetoSols;
    REMOVE.LO = REMOVE.UP;
    REMOVE.L = REMOVE.LO;

    Solve WaterQualityRem USING LP Minimizing obj;

**** record stuff about the solution
    VERTS(v,i,w,s)$((ord(v) eq cV) and CAPP(i,s)) = Areaused.L(i,w,s);
    ObjVal(v,"obj")$(ord(v) eq cV) = Obj.L;
    ObjVal(v,"modstat")$(ord(v) eq cV) = WaterQualityRem.Modelstat;
    ObjVal(v,"solstat")$(ord(v) eq cV) = WaterQualityRem.Solvestat;
    ObjVal(v,"vgen")$(ord(v) eq cV) = cV;
    ObjVal(v,"TotRem")$(ord(v) eq cV) = REMOVE.L;
    );


DISPLAY VERTS, ObjVal, AllDists,vNONZEROVAR, BINDS, HSJstart, NumHSJSols, VertsToUse, cV, lV, NumSolvs;

Parameter Premove(v,w,s,i) Phosphorus removal (kg),
          Pspace(v,w,s,i) Number of spaces need to print phosphorus removal
          PDummy(w,s,i)   Dummy to get listing of watersheds sources and indexes correct
         ;

*Convert area to phosophorus removal
Premove(v,w,s,i) = ef(i)*Verts(v,i,w,s);
Pspace(v,w,s,i) = 3;
Pspace(v,w,s,i)$(Premove(v,w,s,i)) = floor(4+log(Premove(v,w,s,i))/log(10));
PDummy(w,s,i) = 1;

*Write out all the results to GDX

Execute_Unload "WQNE_outG6.gdx";





