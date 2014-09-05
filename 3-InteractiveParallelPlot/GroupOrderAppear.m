function [vUniques, ix, uIndex, iCBeg, iVBeg] = GroupOrderAppear(List)
%   %takes a list of items, groups duplicates, assigns a number to each
%   group and returns the group to which each item belongs as well as the
%   order in which the groups appear in the list.

%   OUTPUTS
%    vUniques = list of groups
%    ix = the last row in List in which each group appears
%    uIndex = the group to which the list member belongs
%    iCBeg = list who order indicates the order the group first appeared in the list 
%    iVBeg = group order associated with the list item

    [vUniques,ix,uIndex] = unique(List);
    lNumGroups = length(vUniques);
    lNumItems = length(List);
    
    blAsCell = iscell(List);
    
     %size and determine where the first member of each group occurs
    if blAsCell
        iCBeg = cell(lNumGroups,1);
        iCBeg{1} = vUniques{uIndex(1)};
    else
        iCBeg = zeros(lNumGroups,1);
        iCBeg(1) = vUniques(uIndex(1));
    end
    
    iVBeg = zeros(lNumItems,1);
    
    cGroup = 1;
    iVBeg(1) = cGroup;
    
    for j=2:lNumItems
        
        [j cGroup];
        if (uIndex(j) == uIndex(j-1))
            iVBeg(j) = iVBeg(j-1);
        elseif (sum((uIndex(j)==uIndex(1:j-1)))==0)
            cGroup = cGroup+1;
            [j cGroup];
            if blAsCell
                iCBeg{cGroup} = vUniques{uIndex(j)};
            else
                iCBeg(cGroup) = vUniques(uIndex(j));           
            end
            iVBeg(j) = cGroup;
        else %fig out which prior group it is
            for k=1:j-1
                if uIndex(j)==uIndex(k)
                    [j k cGroup];
                    iVBeg(j) = iVBeg(k);
                    break;
                end
            end
        end
    end    
end

