function [NewRamp, hWind] = ExpandColorRamp(n,OrigRamp,OrigSpace,ShowGraph)
%   Expands the original color ramp to include n intermediate values
%   between 1st and OrigSpace+1 entries in the ramp. Uses linear
%   interpolation on each RGB hue. Will repeat for the OrigSpace+2 and
%   2*OriginalSpace+1 entries, etc.
%
%   Orig Ramp is a d by 3 matrix of d colors and R G B in the 3 succesive
%   columns.
%   OrigSpace is the original spacing of start -> end elements in the
%   OrigRamp.
%   If ShowGraph is true or positive, then show the graph as well
%   Set OrigSpace = 0 to rather linearly interpret n values by nearest neighbor

%   Returns the NewRamp of colors (d+n*g by 3) as well as the handle hWind
%   to the bar graph created

    [d c] = size(OrigRamp);

    if c~=3
        error('ExpandColorRamp: OrigRamp must have 3 columns.'); 
    end
    
    if OrigSpace == 0
        %Interpret using nearest neighbor
                
        vNs = [1:(d-1)/(n-1):d];
        lNs = length(vNs);
        NewRamp = zeros(lNs,3);
        
        for g=1:n
            startP = floor(vNs(g));
            endP = ceil(vNs(g));
            
            if startP==endP
                NewRamp(g,:) = OrigRamp(startP,:);
            else
                for k=1:3
                    NewRamp(g,k) = interp1([startP endP],OrigRamp([startP endP],k),vNs(g),'linear');
                end
            end
        end
        mPlot = ones(2,lNs);
    else

        nG = floor(d/(OrigSpace+1));
        %Interpolate the ramp

        NewRamp = zeros((n+2)*(nG),3);

        for g=1:nG
            i=(OrigSpace+1)*(g-1)+1;
            j=(n+2)*(g-1)+1;

            [g i j];

            for k=1:3
                if OrigRamp(i,k) == OrigRamp(i+OrigSpace,k)
                    NewRamp(j:j+n+1,k) = OrigRamp(i,k);
                else
                    NewRamp(j:j+n+1,k) = [OrigRamp(i,k):(OrigRamp(i+OrigSpace,k)-OrigRamp(i,k))/(n+1):OrigRamp(i+OrigSpace,k)]';
                end
            end
        end
        
        mPlot = ones(2,(n+2)*nG);
    end
    
    %Plot up the results as a stacked bar graph
    if ShowGraph>0
        Fig1 = figure;
        bh = bar([1:2],mPlot,'Stacked');

        for i=1:length(bh)
            set(bh(i),'facecolor',NewRamp(i,:));
        end

        hWind = bh;
    else
        hWind = 0;
    end
end

