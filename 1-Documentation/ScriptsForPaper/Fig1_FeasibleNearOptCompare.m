function [Fig1,Fig2] = Fig1_FeasibleNearOptCompare(c, A, rel, b, direction, tolerance, fontsize, axislabels, contourspace, omitNE, MGAInfo, blManualContour)
% Graphs the feasible region and near-optimal regions
% for the Convex Nonlinear otimization program
% min (max)z = f(X) = c1*X1^2 + c2*X1 + c3*X2^2 + c4*X2 + c5
% % Subject to Ax <= b (or Ax >= b),
% x >= 0
%
% Problem is small and the near-optimal region is identified by vertex
% enumeration.
%
% This is the script used to generate Figure 1 in the paper.
% EXAMPLE USE:
% [fig1, fig2] = Fig1_FeasibleNearOptCompare([1 -20 10 0 100], [3 2;0 1], ['<' '<'], [12 4]','min',1.8,14, ...
%         {'Decision Variable 1 (X_1)' 'Decision Variable 2' '(X_2)'},20,0,[1 0 15;2 0 1],0);


% INPUTS
%   rel is a cell array that describes the direction of the constraints '<' or '>'
%   direction = 'max' or 'min' and describes whether the problem is a
%        minimization or maximimization problem.
%
%   tolerance = expresses the near optimal region and a fraction (< 1 for a maximization program; > 1 for a minimization problem) of the
%      optimal objective function value
% 
%      near optimal adds a series of piecewise linear constraints
%      c1x1 + x2x2 >= f*(tolerance) ((for a maximimization problem))

%   fontsize = the size to print primary (axis, title, etc. fonts)op

%   axislabels = 2 element cell array for the labels of Axis 1 and Axis 2

%   contourspace = spacing of the objective function contours

%   omitNE = boolean variable that takes the value of 1 to omit the near optimal region and
%       shows only the feasible region. If 0 or omited, shows the Near Optimal
%       region
%
%   MGAInfo = an i x 3 matrix that specifies the i method(s) to use to
%       generate MGA alternatives to include on the plot and associated
%       parameters for each method. Default value is a single scalar 0 (don't
%       include MGA methods).
%       
%          Column 1 - the MGA method to use(See further descriptions in
%                   doMGA.m)
%               0 = ignore this row
%               1 = By Hop-Skip-Jump
%               2 = Maximally different by serial
%               3 = Maximally different all at once
%           Column 2 - Solve by: enter 1 to solve with the Matlab genetic algorithm, 0 = with
%               classic LP or NLP optimization (Default: 0)
%           Columns 3 - maximum number of alternatives to generate
%               (Default: 10)
%
%   blManualContour = boolean takes the value of 1 to let the user manually
%       set contour labels on the plot. 0 or omitted, automatically set
%       contour labels.
%
% OUTPUTS
%   Fig1 = handle to the first plot generated (cartesian plot of the example problem)
%   Fig2 = handle to the 2nd plot generated (parallel coordinate plot of
%           the example problem)
%
% CALLED FUNCTIONS
%    - toolbox:  \optim\optim\linprog.m
%    - other:  delcols
%    - other:  extrdir
%    - other:  extrpts
%    - other:  polygeom
%    - other:  doMGA (generate MGA alternatives)
%    - other:  stratgibbs (generate near-optimal alternatives by stratified
%                   sampling)
%    - other:  nearoptplotmo2 (plot the results in parrallel coordinates)
 
% 
% #####################
%   Programmed by David E. Rosenberg
%   July, 2011 (for career proposal)
%   Sep, 2013 -- updated to fix for paper
%   July, 2014 -- updated to include Modelleing to Generate Alternative
%   points on the plot as defined by Brill et al. 1982. 
%   December, 2014 -- updated to allow the user to specify the type of MGA
%       alternatives to add to the plot. Also, show the parallel cooridnate
%       plot for the same data.
%      
%   Dept. of Civil & Env. Engineering and Utah Water Research Lab
%   Utah State University
%   david.rosenberg@usu.edu
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
%   Bug Reports:
%   Bug reports are much appreciated. Please submit via the the issue tracker on the
%   GitHub repository where you downloaded this file.
%   Note, while much appreciated, there is no promise of when, if the bug will be corrected.

    %%
    % Error checking on user provided inputs
    [m, n] = size(A);
    if n ~= 2
        str = 'Number of the legitimate variables must be equal to 2';
        msgbox(str,'Error Window','error')
        return
    end


    if size(c,2) ~= 5
        str = 'Number of the objective function coefficients must be equal to 5';
        msgbox(str,'Error Window','error')
        return
    end

    %First do the feasible region
    [vert bases] = extrpts(A,rel,b);

    CornerPts = vert';

    if isempty(vert)
        disp(sprintf('\n Empty feasible region'))
        return
    end
    vert = vert(1:2,:);
    vert = delcols(vert);
    d = extrdir(A,b,rel);
    if ~isempty(d)
        msgbox('Unbounded feasible region','Warning Window','warn')
        disp(sprintf('\n Extreme direction(s) of the constraint set'))
        d
        disp(sprintf('\n Extreme points of the constraint set'))
        vert
        return
    end

    if nargin<11
        MGAInfo = 0;
        mgaSz=0;
    else
        %Screen out rows with a method of zero   
        MGAInfo = MGAInfo(MGAInfo(:,1)>0,:);
        [mgaSz,mgaCol] = size(MGAInfo);

        if mgaSz ==0
            MGAInfo = 0;
            warning('No MGA methods entered. Continuing without.')
        else
            if mgaCol < 2
                %add default values for columns 2 and 3
                MGAInfo = [MGAInfo zeros(mgaSz,1) 10*ones(mgaSz,1)];
            elseif mgaCol < 3
                %add default values for column 3
                MGAInfo = [MGAInfo 10*ones(mgaSz,1)];
            end

            hMga = zeros(1,mgaSz); %Handles for the plot markers for alternatives generated by each method
        end
    end

    feColor = [1 .316 1];

    t1 = vert(1,:);
    t2 = vert(2,:);
    %test = 'Before Finished ConvHull'
    z = convhull(t1,t2);
    %test = 'Finished ConvHull'
    Fig1 = figure;
    hold on
    h1 = patch(t1(z),t2(z),feColor);
    set(h1,'LineWidth',2)
    h = 0;
    plotextra = 0.5;
    mit1 = min(t1)-h;
    mat1 = max(t1);
    mit2 = min(t2)-h;
    mat2 = max(t2);

    %find the optimal point

    zvalues = c(1)*(t1(:)).^2 + c(2)*t1(:) + c(3)*(t2(:)).^2 + c(4)*t2(:) + c(5);
    %set the header to show the optimization problem
    if direction=='max'
        str = 'Maximize ';
        zopt = max(zvalues);
    else
        str = 'Minimize ';
        zopt = min(zvalues);
    end

    nopt = tolerance*zopt;

    k=max(size(zvalues));
    %find the optimal decisions
       for j=1:k
            if zopt==zvalues(j)
                x1opt=t1(j);
                x2opt=t2(j);
            end
        end

    %set up the equalities to display
    equalities = cell(m,1);
    for i=1:m
       if rel(i) == '>'
           equalities(i,1)={'\geq'};
       else
           equalities(i,1)={'\leq'};
       end
    end

    %figure out how to print objective function
    %first term

    appends = {'(X_1)^2' 'X_1' '(X_2)^2' 'X_2' ''};
    finaltxt = cell(5);

    for i=1:5
        if c(i)==0
            finaltxt{i} = '';
        else
            if c(i)==1
                stemp = '+';
            else
                stemp=sprintf('%+0.f',c(i));
            end
            finaltxt{i} = sprintf('%s%s',stemp,appends{i});
        end
    end

    objtext = sprintf('f(X_1,X_2) = %s%s%s%s%s\n%s',finaltxt{1},finaltxt{2},finaltxt{3},finaltxt{4},finaltxt{5},'such that:  ');
    for i=1:m
        if A(i,1) == 0
            s1='';
        elseif A(i,1) == 1
            s1='X_1';
        else
            s1=sprintf('%.0fX_1',A(i,1));
        end
        if A(i,2) == 0
            s2='';
        elseif A(i,2) == 1
            s2='X_2';
        else
            s2=sprintf('%.0fX_2',A(i,2));
        end   


        if size(s2,2)==0 || size(s1,2)==0
            s3='';
        else
            s3='+';
        end

        objtext = sprintf('%s%s %s %s %s %.0f',objtext,s1,s3,s2,char(equalities(i,1)),b(i));
        if i < m
            objtext = sprintf('%s\n',objtext);
        end
    end


    objtext = sprintf('%s%s',objtext,', and X_1,X_2 \geq 0');
    title(sprintf('%s%s',str,objtext))
    axis([mit1 max(mat1,mat2)+plotextra mit2 max(mat1,mat2)+plotextra+1])
    h = get(gca,'Title');
    set(h,'FontSize',fontsize)

    if (nargin >= 9) && all(~strcmpi('',axislabels))
        xtext = axislabels{1};
        ytext = sprintf('%s\n%s',axislabels{2},axislabels{3});
        ytextpar = sprintf('%s %s',axislabels{2},axislabels{3});
    else
        xtext = sprintf('%s\n%s','Decision','Variable 1 (X_1)');
        ytext = sprintf('%s\n%s','Decision','Variable 2 (X_2)');
        ytextpar = sprintf('%s\n%s','Decision','Variable 2 (X_2)');
    end
    xlabel(xtext)
    h = get(gca,'xlabel');
    set(h,'FontSize',fontsize)
    ylabel(ytext)
    h = get(gca,'ylabel');
    set(h,'FontSize',fontsize)
    set(gca,'box','on');
    grid

    currFig = get(h,'Parent');

    IncludeNE = (nargin<10) || ((nargin>=10) && (omitNE==0));

    %label the feasible region
    [geom, iner, cpmo] = polygeom(t1,t2); % Find the centroid of the region to center the label on
    textFesReg = text(geom(2), geom(3),'Feasible region','FontSize',fontsize,  ...
           'FontWeight', 'normal', 'HorizontalAlignment', 'center');

    %near optimal region
    if (tolerance <= 0)
        disp(sprintf('Tolerance = %.2f and must be greater than 0',tolerance))
        return
    end 

    if strcmp(direction,'min') && (tolerance < 1)
        disp(sprintf('Tolerance = %.2f and must be greater than 1 for a minimization problem',tolerance))
        return    
    end

    if strcmp(direction,'max') && (tolerance > 1)
        disp(sprintf('Tolerance = %.2f and must be less than 1 for a maximization problem',tolerance))
        return    
    end

    interppoints = 10;

    %choose a dimension, X1 and find the range over which the variable can vary
    x1range = [min(t1) max(t1)];
    x1int = [x1range(1):(x1range(2)-x1range(1))/interppoints:1.7 1.8:0.05:2.1 2.3:(x1range(2)-x1range(1))/interppoints:x1range(2)];
    %x1int = [x1range(1):(x1range(2)-x1range(1))/interppoints:x1range(2)];
    %find the corresponding points in x2 by reworking the objective function

    [j k] = size(x1int);
    x2int=zeros(j,k);
    complex=zeros(j,k);

    for i=1:k
        quad_a = c(3);
        quad_b = c(4);
        quad_c = -(nopt - c(5) - c(1)*x1int(i).^2 - c(2)*x1int(i));

        if quad_a == 0
            x2int(i) = -quad_c/quad_b;
        elseif quad_b^2 - 4*quad_a*quad_c > 0
            x2int(i) = (-quad_b + (quad_b^2 - 4*quad_a*quad_c)^0.5)/(2*quad_a);

        %if nopt > c(1)*(x1int(i)-c(3))^2
        %    x2int(i) = (((nopt - c(1)*(x1int(i)-c(3))^2))/c(2))^0.5 + c(4);

            %check if veritcal increase is too much
        else
           complex(i)=1;
        end
    end

    %add the points to the k-1 pts to the matrix; each represents a line
    %between adjacient points
    %[x1int' x2int']
    Anearopt = A;
    bnearopt = b;
    relno = rel;

    for i=1:k-1
       if complex(i)==0
           m = (x2int(i+1)-x2int(i))/(x1int(i+1)-x1int(i));
           Anearopt = [Anearopt; -m 1];
           bnearopt = [bnearopt; x2int(i) - m*x1int(i)];
           if direction=='max'
               relno = [relno '>'];
           else
               relno = [relno '<'];
           end
       end
    end

    [vert bases] = extrpts(Anearopt,relno,bnearopt); 

    if isempty(vert)
        disp(sprintf('\n Empty feasible region'))
        return
    end

    vert = vert(1:2,:);
    vert = delcols(vert);
    d = extrdir(Anearopt,bnearopt,relno);
    if ~isempty(d)
        msgbox('Unbounded feasible region','Warning Window','warn')
        disp(sprintf('\n Extreme direction(s) of the constraint set'))
        d
        disp(sprintf('\n Extreme points of the constraint set'))
        vert
        return
    end
    t1 = vert(1,:);
    t2 = vert(2,:);

    z = convhull(t1,t2);
    hold on

    neColor = [.316 1 .316]; %pink for near-optimal color

    sOptText = sprintf('f*(%.0f,%.0f)=%.0f',x1opt,x2opt,zopt);
    
    vDistStore = {};
    lDistLogged = [];

    if IncludeNE

        hNearOpt = patch(t1(z),t2(z),neColor);
        set(hNearOpt,'LineWidth',2)
        strfeas = sprintf('%s\n%s','Remaining', 'feasible region');
        set(textFesReg,'string',strfeas);

        h = 0;

        vResults = []; %for comparing results across methods

       if mgaSz > 0
          %Calculate the MGA alternatives
          xMGA = [x1opt x2opt]; %seed the first solution as the optimum
          %Add non-negativity constraints
          AnearoptMGA = [Anearopt;-1 0; 0 -1];
          bnearoptMGA = [bnearopt; 0;0];

          xOpt = xMGA;
          %Nonlinear programming formulation
          lB = zeros(2,1);
          uB = [4;5];

          sMarkers = {'s' '^' '+' 'x' '*'};
          mgaNames = {'MGA-HSJ' 'MGA-Serial' 'MGA-Simultaneous'};
          solveBy = {'Opt' 'GA'};
          vMarkerSize = [10 10 9 10 10];

          sColor = zeros(5,3); %[0 0 0;0.316 0 0.316; 0.526 0 0.526];
          sColor(3,:) = [0.6 0.06 0.06];

          %Loop through the alternative generation methods
          for i=1:mgaSz

             mgaNames{MGAInfo(i,1)}

             %Linear programming formulation
             [xMGA, p, MinDist, vDists, RetFlag, RetStr, rTime] = doMGA(xOpt,struct('A',AnearoptMGA,'b',bnearoptMGA,'MaxAlts',MGAInfo(i,3),'StopTime',60,'MGAType',MGAInfo(i,1),'errorcrit',1e-6,'StopDistance',0.2,'SolveAsGA',MGAInfo(i,2)));
             %Non-linear programming formulation
             %[xMGA, p, vDists, RetFlag, rTime] = doMGA(xOpt,struct('A',A,'b',b,'lB',lB,'uB',uB,'MaxAlts',10,'MGAType',i,'StopDistance',0.2, ...
             %                         'objfunc',@(x)c(1)*(x(1)).^2 + c(2)*x(1) + c(3)*(x(2)).^2 + c(4)*x(2) + c(5),...
             %                          'tolerance',tolerance));

             vResults = [vResults; mgaNames{MGAInfo(i,1)} solveBy{MGAInfo(i,2)+1} num2cell(p) num2cell(MinDist) num2cell(rTime) num2cell(RetFlag) RetStr];

             if ~isempty(vDists) && (p>1)
                lDistLogged = [lDistLogged; i];
                vDistStore = {vDistStore{:} vDists};
             end

             uistack(Fig1); %Bring to the front for further plotting

              if size(xMGA,1)>0

                  %Plot the MGA point(s) as purple triangles
                  hMga(i) = plot(xMGA(:,1),xMGA(:,2));
                  set(hMga(i),'linestyle','none','marker',sMarkers{i},'markerfacecolor',sColor(i,:),'color',sColor(i,:),'markersize',vMarkerSize(i));
                  %text(xMGA(2,1),xMGA(2,2),sprintf(mgaNames{i}), 'HorizontalAlignment', 'right','VerticalAlignment','bottom','color',[0 0 0],'FontSize',fontsize-2)
                  if MGAInfo(i,1) == 0 %ignore
                      %Label the MGA points in the order they were generated
                      for j=1:size(xMGA,1)
                          text(xMGA(j,1),xMGA(j,2),sprintf('%d',j),'color',[0 0 0],'FontSize',fontsize-2)
                      end
                  end
              end
          end

          %Stack the first one back on top
          uistack(hMga(1));

       else    

        %Find the centroid of the polygon so we can label it
        [geom, iner, cpmo] = polygeom(t1,t2);
        text(geom(2), geom(3),sprintf('%s\n%s\n%s','Near-', 'optimal','region'),'FontSize',fontsize,  ...
           'FontWeight', 'normal', 'HorizontalAlignment', 'center')

        sOptText = sprintf('Optimum:\n%s',sOptText);
       end
    end

    %Plot the optimum point
    %plot optimal point as blue circle
    hOpt = plot(x1opt,x2opt);
    optColor = [0 0 1];
    set(hOpt,'linestyle','none','marker','o','markerfacecolor',optColor,'color',optColor,'markersize',10);

    text(x1opt,x2opt,sOptText, 'HorizontalAlignment', 'left','VerticalAlignment','bottom','color',optColor,'FontSize',fontsize-2)

    %plot objective contours
    %generate the mesh
    mmin = min(mit1,mit2);
    mmax = max(mat1+plotextra,mat2+plotextra+1);
    x1 = [mmin:(mmax-mmin)/contourspace:mmax];

    [X1,X2] = meshgrid(x1);

    Z = c(1)*(X1).^2 + c(2)*X1 + c(3)*(X2).^2 + c(4)*X2 + c(5);
    %Z = c(1)*(X1-c(3)).^2+c(2)*(X2-c(4)).^2;

    zRange = [min(min(Z)) max(max(Z))];

    grid off

    cLabels = {'Optimum' 'Near-optimal'};

    if mgaSz > 0
        %include a legend
        hLegend = legend([hOpt; hMga(hMga>0)'; hNearOpt],{cLabels{1} mgaNames{MGAInfo(hMga>0,1)} cLabels{2}},'FontSize',fontsize-2);
    end


    [h1,cs] = contour(X1,X2,Z,[0:contourspace:1000]);
    %Label the contours
    if (nargin>11) && (blManualContour == 1)
        clabel(h1,cs,'manual'); %'LabelSpacing',72*50,
    else
        clabel(h1,cs,'LabelSpacing',72*3)
    end

    set(cs,'color',optColor,'LineStyle','-.');

    hold off


    % Plot the parallel coordinate version of the objective function and
    % feasible region in a second plot. Optimal solution in row 1, random
    % sampled alternatives in subsequent rows

    vDirection = 2*strcmpi(cellstr(rel'),'<')-1; %map 1/0 less than/greater than to 1/-1

    Asamp = repmat(vDirection,1,size(A,2)).*A;
    bsamp = vDirection.*b;

    %Append non-negativity constraints
    Asamp = [Asamp; -eye(2)];
    bsamp = [bsamp; zeros(2,1)];

    NewP = [CornerPts; vert'];
    vValid = ones(size(NewP,1),1);   
    
    %Build the optimization model formulation to pass to the sampling
    %routine
      ModelForm.Aineq = Asamp;
      ModelForm.bineq = bsamp;
      %ModelForm.lb = zeros(2,1); %Non-negativity constraints
      ModelForm.solver = 'linprog';
      ModelForm.options = struct('Algorithm','simplex','Display','off');
    

    [SampledP,vSampValid, sgRunTime] = stratgibbs(2500,ModelForm,struct('extmethod','opt','x0',[x1opt x2opt],'errorresid',0,'GibbsDrawsPerSample',1));

    NewP = [NewP; SampledP];
    vValid = [vValid; vSampValid];

    %hold on
    %plot(NewP(:,1),NewP(:,2),'line','none','marker','+');
    Fig2 = [];

    if 1 % exclude 2nd plot
        hold off

        NewZ = c(1)*(NewP(vValid>0,1)).^2 + c(2)*NewP(vValid>0,1) + c(3)*(NewP(vValid>0,2)).^2 + c(4)*NewP(vValid>0,2) + c(5);

        mData = [zopt x1opt x2opt; NewZ NewP(vValid>0,:)];
        vGroups = [cLabels(1); repmat({'Feasible region'},size(mData,1)-1,1)];

        %Identify near-optimal latenratives as a new group

        vGroups(mData(:,1) <= nopt) = cLabels(2);
        vGroups(1) = cLabels(1);

        %Build the colors for the groups so they match figure 1
        mColors = zeros(3,2,3,3);
        mColors(1,1,:,:) = repmat(feColor,3,1);
        mColors(2,1,:,:) = repmat(neColor,3,1);
        mColors(3,1,:,:) = repmat(optColor,3,1);
        mColors(:,2,:,:) = mColors(:,1,:,:);

        %Define the optional parameters needed to plot in parallel coordinates
        vParams = {'FontSize',20,'GroupToHighlight',2,'vGroup',vGroups,'mGroupData',[{'Feasible region'} 1 1; {'Near-optimal'} 1 1; {'Optimum'} 1 1.5], ...
               'vObjLabels',{'Cost'},'vXLabels',{xtext ytextpar},'yAxisLabels',{'Objective function value' 'Decision variable value'},...
               'AxisScales','custom',[1 1],[0 0 0; 300 5 5],'NumTicks',5, ...
               'StartTab',3,'ShowObjsDiffColor',0,'ShowGroupLabels',0,'HideCheckboxes',1 ...
               'ShowControls',0,'ShowInsetPlot',0,'GroupToHighlight',2,'mHighlightColor',optColor,...
               'mColors',mColors,'PlotPosition',[0.100 0.3 0.47 0.6],'mColorYScales',[0 0 0;0 0 0],'NumTicks',[6 6],'YAxisMargins',[0.1 0.1]};
        Fig2 = nearoptplotmo2(mData(:,1),mData(:,2:end),vParams{:});

        %Append headers plus results for random sample to compare across near-opt methods

        vResults = [vResults; 'Random sample' ' ' num2cell(sum(vSampValid>0)) num2cell(0) num2cell(sgRunTime) num2cell(1) sprintf('%d of %d valid',size(vSampValid>0,1),size(vSampValid,1))];
    end

    %Print out summary of results
    for i=1:length(lDistLogged)
        mgaNames{MGAInfo(lDistLogged(i),1)}
        {'Alt' 'Min. Distance' 'Cum. Run Time (sec)'}
        [num2cell([1:size(vDistStore{i},1)]') num2cell(vDistStore{i})]
    end

    [{'Gen. method' 'By' '# Alts' 'Min Dist' 'Run Time (sec)' 'Flag' 'Explanation'}; vResults]


end
