function [c,r] = chebycenterFull(probstruct)
%CHEBYCENTER Compute Chebyshev center of polytope with the problem structure
%  defined by problstruct with inequality, equality, lower-bound, and upper-
%  bound constraints as well as options.
%  The Chebyshev center of a polytope is the center of the largest
%  hypersphere enclosed by the polytope. 
%  Requires optimization toolbox.
%  Allows equality, lower-bound and upper-bound constraints -- DER February 2015
%  Keep lower- and upper-bound constraints in original form because they
%  solve cleaner

% Check the problem structure
    [a_full,b_full] = OptimiFull(probstruct);
    
    if isempty(a_full)
        warning('Problem with problem structure')
        return
    end

%Defaults
A1ineq = []; bineq=[];
A1eq = []; beq = [];
lb1 = []; ub1 = [];

if isfield(probstruct,'Aineq') && isfield(probstruct,'bineq')
    Aineq = probstruct.Aineq;
    bineq = probstruct.bineq;
    
    [n,p] = size(Aineq);
    an = sqrt(sum(Aineq.^2,2));
    A1ineq = [Aineq an];
end

if isfield(probstruct,'Aeq') && isfield(probstruct,'beq') && ~isempty(probstruct.Aeq) && ~isempty(probstruct.beq)
    Aeq = probstruct.Aeq;
    beq = probstruct.beq;
    [nEq,p] = size(Aeq);
    %aneq = sqrt(sum(Aeq.^2,2));
    A1eq  = [Aeq zeros(nEq,1)];
end

if isfield(probstruct,'lb') && ~isempty(probstruct.lb)
    lb = probstruct.lb;
    p = size(lb,1);
    lb1 = [lb;sqrt(sum(lb.^2))];
end
if isfield(probstruct,'ub') && ~isempty(probstruct.ub)
    ub = probstruct.ub;
    p = size(ub,1);
    ub1 = [ub;sqrt(sum(ub.^2))];
end

    f = zeros(p+1,1);
    f(p+1) = -1;

[c, fopt, exitflag] = linprog(f,A1ineq,bineq,A1eq,beq,lb1,ub1,[],probstruct.options);
%c = linprog(f,A1,b,[],[],[],[],[],options);

if exitflag ~=1
    %exitflag
    c=[];
    r=NaN;
else
    r = c(p+1);
    c = c(1:p);
end
end
