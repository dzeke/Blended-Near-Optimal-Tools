function Fig1_FeasibleNearOptCompare(c, A, rel, b, direction, tolerance, fontsize, axislabels, contourspace, omitNE, blAddMGA)
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
% Fig1_FeasibleNearOptCompare([1 -20 5 0 100], [3 2;0 1], ['<' '<'], [12 4]','min',1.8,14, {'Decision Variable 1 (X_1)' 'Decision Variable 2' '(X_2)'},20,1,1)

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

%   fontsize = the size to print primary (axis, title, etc. fonts)

%   axislabels = 2 element cell array for the labels of Axis 1 and Axis 2

%   contourspace = spacing of the objective function contours

%   omitNE = boolean variable that takes the value of 1 to omit the near optimal region and
%       shows only the feasible region. If 0 or omited, shows the Near Optimal
%       region
%
%   blAddMGA = boolean takes the value of 1 to calculate and include MGA
%       maximally different solutions on the plot. 0 == don't plot MGA
%       solutions. Default value: 0.
%
% CALLED FUNCTIONS
%    - toolbox:  \optim\optim\linprog.m
%    - other:  delcols
%    - other:  extrdir
%    - other:  extrpts
%    - other:  polygeom
 
% 
% #####################
%   Programmed by David E. Rosenberg
%   July, 2011 (for career proposal)
%   Sep, 2013 -- updated to fix for paper
%   July, 2014 -- updated to include Modelleing to Generate Alternative
%   points on the plot as defined by Brill et al. 1982. 
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
%First do the feasible region
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


[vert bases] = extrpts(A,rel,b);

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
    blAddMGA = 0;
end

t1 = vert(1,:);
t2 = vert(2,:);
%test = 'Before Finished ConvHull'
z = convhull(t1,t2);
%test = 'Finished ConvHull'

hold on
h1 = patch(t1(z),t2(z),[1 .316 1]);
set(h1,'LineWidth',2)
h = 0;
plotextra = 0.5;
mit1 = min(t1)-h;
mat1 = max(t1);
mit2 = min(t2)-h;
mat2 = max(t2);

%h2 = axes('Position',get(h1,'Position'));

%set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])
%set(h2,'XLim',get(h1,'XLim'),'Layer','top')


%set(h,'linestyle','--','color',[0 0 0])
%set(h,'linewidth',2)

%find the optimal point

%zvalues = c(1)*(t1(:)-c(3)).^2+c(2)*(t2(:)-c(4)).^2;
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

if nargin >= 9
    xtext = axislabels{1};
    ytext = sprintf('%s\n%s',axislabels{2},axislabels{3});
else
    xtext = 'Decision variable 1 (X_1)';
    ytext = sprintf('%s\n%s','Decision variable 2','(X_2)');
end
xlabel(xtext)
h = get(gca,'xlabel');
set(h,'FontSize',fontsize)
ylabel(ytext)
h = get(gca,'ylabel');
set(h,'FontSize',fontsize)
grid

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
[x1int' x2int']
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

%Anearopt
%bnearopt
%relno

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



if IncludeNE
    
    h = patch(t1(z),t2(z),[.316 1 .316]);
    set(h,'LineWidth',2)
    strfeas = sprintf('%s\n%s','Remaining', 'feasible region');
    set(textFesReg,'string',strfeas);

    h = 0;

    %Find the centroid of the polygon so we can plot it
    [geom, iner, cpmo] = polygeom(t1,t2);
    text(geom(2), geom(3),sprintf('%s\n%s\n%s','Near-', 'optimal','region'),'FontSize',fontsize,  ...
       'FontWeight', 'normal', 'HorizontalAlignment', 'center')
   
   if blAddMGA
      %Calculate the MGA alternatives
      xMGA = [x1opt x2opt]; %seed the first solution as the optimum
      NonZeroDecs = (xMGA>0)+0; %The non-zero deicision variables. Take largest value of any prior solutions. Once a decision is non-zero, always non-zero. 0 at end to turn to numeric
      kMGA = 1; %index of MGA alternatives
      NumNonZerosDecs = sum(NonZeroDecs,2);
      
      %Add non-negativity constraints
      AnearoptMGA = [Anearopt;-1 0; 0 -1];
      bnearoptMGA = [bnearopt; 0;0];
      
      %Greater = strcmpi(relno','>');
      %(Greater==1,:) = -Anearopt(Greater==1,:);
      %bnearoptMGA(Greater==1) = -bnearopt(Greater==1);
      
      while (NumNonZerosDecs(kMGA)<n) || ((kMGA>1) && (NumNonZerosDecs(kMGA) > NumNonZerosDecs(kMGA-1)))
           
          [xNew, fNew, errorflag, Test2] = linprog(NonZeroDecs(kMGA,:),AnearoptMGA,bnearoptMGA);
           
          if errorflag==1
            %Record the solution
            xMGA = [xMGA;xNew'];
 
            kMGA = kMGA+1;
            
            NonZeroDecs(kMGA,:) = (max(xMGA)+0>0); %
            NumNonZerosDecs(kMGA) = sum(NonZeroDecs(kMGA,:),2);
          else
              Test2.message
              warning('MGA terminiated earlier')
          end
      end 
      %Plot the MGA point(s) as purple triangles
      hMga = plot(xMGA(2:end,1),xMGA(2:end,2));
      set(hMga,'linestyle','none','marker','^','markerfacecolor',[0.316 0 0.316],'color',[0.316 0 0.316],'markersize',8);
      text(xMGA(2,1),xMGA(2,2),sprintf('MGA'), 'HorizontalAlignment', 'right','VerticalAlignment','bottom','color',[0.316 0 0.316],'FontSize',fontsize-2)

   end
end

%Plot the optimum point
%plot optimal point as black circle
h = plot(x1opt,x2opt);
%set(h,'linestyle','none','marker','o','markerfacecolor',[0 0.526 0],'color',[0 0.526 0]);
set(h,'linestyle','none','marker','o','markerfacecolor',[0 0 1],'color',[0 0 1],'markersize',8);
text(x1opt,x2opt,sprintf('Optimum:\nf*(%.0f,%.0f)=%.0f',x1opt,x2opt,zopt), 'HorizontalAlignment', 'left','VerticalAlignment','bottom','color',[0 0 1],'FontSize',fontsize-2)

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

[h1,cs] = contour(X1,X2,Z,[0:contourspace:1000]);
%Label the contours - instead use the commented line #428 to manually
%   position labels on the plot
%clabel(h1,cs,'manual'); %'LabelSpacing',72*50,
clabel(h1,cs,'LabelSpacing',72*3)
set(cs,'color',[0 0 1],'LineStyle','-.')

hold off
