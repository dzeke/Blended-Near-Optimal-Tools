function [print_form] = PrintGamsSo(inp_so,InputStruct)
% Formats a GAMS structure in cell/matrix form so that UELs display as text to indicate the value

% INPUTS
% inp_so = the input structure. Must have val and uel fields
% InputStruct = optional paramter takesthe value of 1 to indicate inp_so is
%               a structure passed to GAMS (with separately indexed.val and .uels fields)
%               Default value is 0 meaning inp_so is a structure passed
%               from GAMS back to Matlab (with .val and .uels fields with
%               the first 1 to n-1 columns indexes into the uels.
% 
% OUTPUTS
% print_form = final output matrix

    if (nargin > 1) && (InputStruct==1)
        sz_inp = size(inp_so.uels);
        if sz_inp(2) > sz_inp(1)
            print_form = [inp_so.uels' num2cell(inp_so.val')];
        else
            print_form = [inp_so.uels num2cell(inp_so.val)];
        end
    else

         [m n] = size(inp_so.val);
         print_form = [];

         for i=1:n-1
            print_form = [print_form inp_so.uels{1,i}(inp_so.val(:,i))'];
         end

         %last column is the values
         print_form = [print_form num2cell(inp_so.val(:,n))];
    end
end

