function [varargout] = ReorderCols(vNewOrder,varargin)
    %Reorders the columns of the matrixes provided in varargin according to
    %the row vector of column indexes listed in vNewOrder
    
    %Use to reorder or select a subset of columns
    
    %Returns a variable argument list of the same size as varargin and in the same listing order with the
    %reorderd matrices
    
    %All input matrices must have the same number of columns!
    
    %EXAMPLES
    %The examples use the following matricies
    % A = [ones(3,1) 2*ones(3,1) 3*ones(3,1)]
    % B = [3*ones(7,1) 4*ones(7,1) 5*ones(7,1)]
    
    %1. Swap the first and third columns
    %[Anew,Bnew]=ReorderCols([3 2 1],A,B)
    %
    %  
    
    
    %2. Retain the 2nd and 3rd columns but show the 3rd column first
    %[Anew,Bnew]=ReorderCols([3 2],A,B)

    %David E. Rosenberg
    %Utah State University
    %david.rosenberg@usu.edu
    %June 2013
    
    %%
    
    %Check the inputs
    if nargin==1
        error('ReorderCols:ImproperInput','Need to provide a matrix to work on');
    end
    
    n = size(varargin{1},2);
    
    NumMatricies = nargin-1;
    
    %check the inputs all have the same number of columns
    for i=1:NumMatricies
        if (size(varargin{i},2) ~= n) && (size(varargin{i},2) ~= 0)
            error('ReorderCols:ImproperInput','Input Matrix #%d is of different size than first. Reenter.',i)
        end
    end
    %Check the elements of vNewOrder are within the size
    if max(vNewOrder)>n
        error('ReorderCols:ImproperInput','The value %d in vNewRows is larger than the max number of columns (%d) in the input matrices',max(vNewOrder), n)
    end
    
    varargout = varargin;
    %Loop through each input matrix and return the selected columns
    for i=1:NumMatricies
        cMat = varargin{i};
        if isempty(cMat)
            varargout{i} = cMat;
        else
            varargout{i} =  cMat(:,vNewOrder);
        end
    end 
end

