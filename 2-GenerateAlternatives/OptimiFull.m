function [A_full, b_full] = OptimiFull(optimiStructure)
    % Writes the constraints defined in the optimiStructure fields 
    % (.Aineq, .bineq, .Aeq, .beq, .lb, .ub) as a big
    % constraint set of the form A_full x <= b_full.
    % where A_full = m x n matrix and b_full = m x 1 vector.
    %
    % Returns a warning message if there is an inconsistency between
    % arguments and sets A_full and b_full to empty;
    
    % i.e., Aineq x <= bineq
    %       Aeq x <= beq
    %     - Aeq x <= -beq
    %     - eye(n) <= -lB
    %       eye(n) <= uB
    %
    % David E. Rosenberg. February 2015.
      
    A_full = [];
    b_full = [];
    
    % Inequality constraints
    if isfield(optimiStructure,'Aineq') && isfield(optimiStructure,'bineq')
        if ~isempty(optimiStructure.Aineq) && ~isempty(optimiStructure.bineq)
            if (size(optimiStructure.Aineq,1)==size(optimiStructure.bineq,1))
                A_full = [A_full;optimiStructure.Aineq];
                b_full = [b_full;optimiStructure.bineq];
            else
                warning('OptimiFull: .Aineq and .bineq inconsistent sizes')
                return;
            end
        end
    end
    
    [m,n] = size(A_full);
    
    % Equity constraints
    if isfield(optimiStructure,'Aeq') && isfield(optimiStructure,'beq')
        if ~isempty(optimiStructure.Aeq) && ~isempty(optimiStructure.beq)
            if (size(optimiStructure.Aeq,1)==size(optimiStructure.beq,1)) && ...
                    (size(optimiStructure.Aineq,2)==size(optimiStructure.Aeq,2))
                % First term--less than direction of the equality; second term -
                % greater than direction of inequality
                A_full = [A_full;optimiStructure.Aeq; -optimiStructure.Aeq];
                b_full = [b_full;optimiStructure.beq; -optimiStructure.beq];          
            else
                warning('OptimiFull: .Aeq and .beq inconsistent size or number of columns different than .Aineq');
                A_full = [];
                b_full = [];
            end
        end
    end
    
    % Lower bound constraints
    if isfield(optimiStructure,'lb') && ~isempty(optimiStructure.lb)
        if   (size(optimiStructure.lb,1)==n) && ...
                    (size(optimiStructure.lb,2)==1)
            % First less than direction; second greater than direction
            A_full = [A_full;-eye(n)];
            b_full = [b_full;-optimiStructure.lb];
        else
            warning('OptimiFull: .lb different number of rows (decision variables) than columns in .Aineq')
            A_full = [];
            b_full = [];
            return;
        end
    end
    
    % Upper bound constraints
    if isfield(optimiStructure,'ub') && ~isempty(optimiStructure.ub)
        if  (size(optimiStructure.ub,1)==n) && ...
                    (size(optimiStructure.ub,2)==1)
            % First less than direction; second greater than direction
            A_full = [A_full;eye(n)];
            b_full = [b_full;optimiStructure.ub];
        else
            warning('OptimiFull: .ub different number of rows (decision variables) than columns in .Aineq')
            A_full = [];
            b_full = [];
            return;
        end
    end
end

