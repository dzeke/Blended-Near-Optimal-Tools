function [mCounts] = ConsecutiveCount(mMatrix)
%Counts the number of conseculative,common instances of a value in the rows
%and returns a matrix with those values. Pads recurrences with zeros
%
%
% e.g., the matrix A 1 i
%                  A 1 ii
%                  A 2 i
%                  B 1 i
%                  B 2 i
%                  B 2 ii
%                  B 2 iii

%would return                                        [ 3 2 1
%                                                      0 0 1
%                                                      0 1 1
%                                                      4 1 1
%                                                      0 3 1
%                                                      0 0 1
%                                                      0 0 1]
%                                               
% 

   sz = size(mMatrix);
   mCounts = zeros(sz(1),sz(2));
   
   for j=1:sz(2) % work the columns
       lFirstRow = 1; 
       for i=2:sz(1)
           switch class(mMatrix(:,j))
               case 'cell'
                   bSame = strcmpi(mMatrix(lFirstRow,j),mMatrix(i,j));
               case 'double'
                   bSame = (mMatrix(i,j) == mMatrix(lFirstRow,j));
           end
           
           if ~bSame
               %we have a different value
               mCounts(lFirstRow,j) = i-lFirstRow;
               lFirstRow = i;
           end
       end
       
       mCounts(lFirstRow,j) = sz(1)-lFirstRow+1;
       
   end
end

