function [mDataPU, vMultipliers, BaseColOut] = ScaleMultipleAxes(mData,BaseCol,BaseMult)
    %Scales the columns of the m by n matrix of data mData comprising n columns of independent (and differently scaled data sets) 
    %so as to use a common set of ticks and plot units to plot all the
    %columns on the same plot. This scaling allows side-by-side comparison
    %of the dynamic range of data in each column.
    %
    %Also returns:
    % mDataPU = transformed data to plot in units of the dTicks
    % vMultipliers = multiplier applied to mDataPU to get back the original data
    % BaseColOut = the index of the column on which the scaling occurs.
    %   This column will also have the "largest" value and should be used
    %   when constructing a set of ticks for the joint axes
    %
    % Two optional parameters:
    %   BaseCol -- specifies which column to use as the base
    %       units (this column is not transformed). Mutlipliers for all other columns will be chosen so the maximum tick for the other column
    %       If not specified, the program identifies the column to minimize the disparatiy between
    %   BaseMult - specifies the multiplier to apply to use
    %       when calculating and displaying the ticks and vMultipliers (in units of power of 10s, e.g., -1, 0, 1, 2, 3...)
    %       Default value is 0 (no multiplier)
    %
    % Algorithm to use:
    %  2. Find largest value in each column, reduce this value to it's
    %       significant figures (i.e., 100 becomes 1; 240 becomes 24)
    %  3. Compute ratio matrix from largest values. Each entry represents the % of the scale that
    %       be used by the data set in the column when the ticks for the
    %       row are used.
    %  4. Find the minimum value in each row of the matrix.
    %  5. Choose the row in #4 that has the largest value. This is the row
    %       to use as the BaseCol
    %  6. Calculate multipliers for each column.
    %  7. Use the multipliers to calculate the plot units for each column
    %  and ticks
    %
    %
    %David E. Rosenberg
    %Utah State University
    %david.rosenberg@usu.edu
    %June 2013

    [m n] = size(mData);
    nTicks = cell(n,1);
    nTickLims = zeros(2,n);
    Magnitude = zeros(1,n);
    SigTicks = zeros(2,n);
    
    %ax2 = axes('visible','off');
    
    %Step 1. Find the max and min values for each column
    nTickLims(1,:) = max(mData);
    nTickLims(2,:) = min(mData);    
    %Loop through each column
    for i=1:n       
        %Step 2 - cont. Reduce the largest native tick value to it's most signficant
        %units
        if nTickLims(1,i)==0
            Magnitude(i)=0;
        else
            Magnitude(i) = floor(log(nTickLims(1,i))/log(10));
        end
        SigTicks(1,i) = roundn(nTickLims(1,i),Magnitude(i)-1)/(10^Magnitude(i));
        if floor(SigTicks(1,i)) ~= SigTicks(1,i)
            SigTicks(1,i) = 10*SigTicks(1,i);
            Magnitude(i) = Magnitude(i)-1;
        end
    end
    %delete(tPlot); delete(ax2);
    
    nTickLims;
    Magnitude;
    SigTicks;
    
    %Step 3. Compute the ratio matrix
    if (nargin==1) || isempty(BaseCol) || (BaseCol<1) || (BaseCol>n) || isnan(BaseCol)
        mRatios = zeros(n);
        for i=1:n %rows
            for j=1:n %cols
                mRatios(i,j) = SigTicks(1,j)/SigTicks(1,i);
                if mRatios(i,j) > 1
                    mRatios(i,j) = mRatios(i,j)*10^(-ceil(log(mRatios(i,j))/log(10)));
                end
                if mRatios(i,j) <= 0.1
                    mRatios(i,j) = mRatios(i,j)*10;
                end
            end
        end
        %Step 4. Find the minimum value in each row
        MinRatios = min(mRatios,[],2);
        %Step 5. Choose the row with the largest minimum ratio
        mRatios;
        [MaxVal BaseColOut] = max(MinRatios);
    else
        BaseColOut = BaseCol;
    end

    SigTicks;
    vMultipliers = Magnitude;
    %dTicks = nTicks{BaseColOut}/(10^vMultipliers(BaseColOut));
    %Step 6. Calculate the multipliers
    if sum(SigTicks(1,BaseColOut)<SigTicks(1,:)) > 0
        %dTicks = 10*dTicks;
        SigTicks(1,BaseColOut) = 10*SigTicks(1,BaseColOut);
        vMultipliers(BaseColOut) = vMultipliers(BaseColOut)-1;
    end
    
    for i=1:n
        if (i~=BaseColOut) && (vMultipliers(i)>0) && SigTicks(1,i)*10<=SigTicks(1,BaseColOut)
            vMultipliers(i) = vMultipliers(i)-1;
        end
    end
                
    vMultipliers;
    %dTicks;
    
    %Adjust and rescale by the BaseMult
    %minMult = min(vMultipliers);
    
    if nargin<=2
        BaseMult = 0;
    end
    
    vMultipliers = vMultipliers+BaseMult; %-minMult;
    %dTicks = dTicks*10^(BaseMult+minMult);
    
    %Step 7. Calculate the plot units
    mDataPU = mData;
    for i=1:n
        mDataPU(:,i) = mData(:,i)/(10^vMultipliers(i));
    end    
end

