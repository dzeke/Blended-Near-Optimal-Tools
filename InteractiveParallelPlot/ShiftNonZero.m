function [vOut, cRMap] = ShiftNonZero(vVector,Dir)
%   Shifts the non-zero elements of vVector in the direction Dir with no
%   wrapping and returns this new vector as vOut.
%   cRMap is the reverve mapping and says where the ith element in vOut came from
%
%   Dir = +1: right for a row vector, down for a column vector
%   Dir = -1: left for a row vector, up for a column vector
%
%   E.g. [vOut, cRMap] = ShiftNonZero([3 4 0 -2 0  1 4 5],-1)
%   returns 
%           Out = [3     4    -2     0     1     4     5     0]
%           cRMap = [1     2     4     3     6     7     8     5]
%           (3rd and 4th elements are swapped and 5th element is moved to
%               the end)


    %Dir;
    %vVector;
    
    WasRowVector = 0;
    
    if isrow(vVector)
        vVector=vVector';
        WasRowVector = 1;
    end
    
    n=length(vVector);
    
    cRMap = [1:n]';
    
    vBinary = vVector ~= 0;
    
    vDiff = diff(vBinary);
    
    if (sum(vBinary)==0) || (all(vDiff==0)) || (Dir==-1 && all(vDiff(2:end)==0)) || (Dir==1 && all(vDiff(1:end-1)==0))
        %All zero, all checked, or only first or last checking and moving further off the end => no moving!
        vOut = vVector;
        if WasRowVector
            vOut = vOut';
            cRMap = cRMap';
        end
        return
    end        
        
    if Dir==1 %move right, need to right shift the difference cals
        vDiff = [NaN;vDiff];
    end

    %vDiff
    
    StartOn = find(vDiff==1);
    EndOn = find(vDiff==-1);
    
    %Clean these up
    if Dir==-1 %Move Left
        if length(StartOn) > length(EndOn) % we finished on a checked axis. Create a dummy off entry at the end that is the last column
            EndOn = [EndOn;n];
        end
        if StartOn(1) > EndOn(1) %We started with axes that are checked that don't need to move. Shorten EndOn
            EndOn = EndOn(2:end);
        end
    else % Move right
        if (length(StartOn) < length(EndOn)) || (StartOn(1) > EndOn(1)) %The first (left most axes) were checked and need to move right. Pad with the index of the right most column
            StartOn = [1;StartOn];
        end
        if length(StartOn) > length(EndOn) % we finished on a checked axis. Create a dummy off entry at the end that is the last column
            StartOn = [StartOn(1:end-1)];
        end
    end
    %StartOn
    %EndOn
    
    vOut=vVector;
    
    for i=1:length(StartOn)
        vOut(StartOn(i):EndOn(i))=circshift(vVector(StartOn(i):EndOn(i)),Dir);
        cRMap(StartOn(i):EndOn(i))=circshift(cRMap(StartOn(i):EndOn(i)),Dir);
    end
    
    if WasRowVector %convert it back to return
        vOut = vOut';
        cRMap = cRMap';
    end
end

