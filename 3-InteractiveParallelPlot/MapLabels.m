function [i_for] = MapLabels(a,b)
% Returns the mapping for where items in list b appear in list a.
% i_for is a vector of indexes representing the forward mapping of a to b, i.e., a(i_for) = b

    %First check that a and b are of the same type
    if strcmpi(class(a),class(b)) == 0
        error('Classes of inputs a (%s) and b (%s) are different',class(a),class(b))
        return
    end
    
    [m,n] = size(b);
    i_for = NaN(m,n);
    i_rev = i_for;
    
    if strcmpi(class(a),'cell')
        for i=1:m
            for j=1:n
                tmp = find(strcmpi(b(i,j),a));
                if ~isempty(tmp)
                    i_for(i,j) = tmp;
                end
            end
        end
    end
    
    if strcmpi(class(a),'double')
        for i=1:m
            for j=1:n
                tmp = find(b(i,j)==a);
                if ~isempty(tmp)
                    i_for(i,j) = tmp;
                end
            end
        end
    end

        
end

