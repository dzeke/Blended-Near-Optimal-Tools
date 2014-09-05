function [varargout] = OSUColorRamps(varargin)
% Returns a n by 3 matrix of colors for the OSU color ramps, one for each
% ramp name provided. The single parameter 'all' for 'varargin' returns all
% the ramps in order as the varargout.

% The color ramps are available at for viewing at
% http://geography.uoregon.edu/datagraphics/color_scales.htm.

% David E. Rosenberg
% Utah State University
% June 2013
% david.rosenberg@usu.edu

%%  Ramp Definitions
    RampNames = {'GreenToMagenta16Step' 'Categorical12Step' 'SteppedSequential' 'BlueToGrey8Step' 'LightToDarkBlue10Step' 'BlueToRed18Step'  'LightToDarkBlue7Step'};

    %GreenToMagenta-16Step
    %   Row 1 is Dark Green, Row 16 is Magenta, Rows 8 and 9 and near white
    GreenToMagenta16Step = [0	0.316	0  %Dark Green
        0	0.526	0
        0	0.737	0
        0	0.947	0
        0.316	1	0.316
        0.526	1	0.526
        0.737	1	0.737
        1	1	1               % Near White
        1	0.947	1           % Near White
        1	0.737	1
        1	0.526	1
        1	0.316	1
        0.947	0	0.947
        0.737	0	0.737
        0.526	0	0.526
        0.316	0	0.316];  % Magenta

    %Categorical 12 Step
    %  Odd rows are light, even rows are dark. Rows 1 and 2 are light and dark Orange
    % Then Yellow, Green, Blue, Purple, and Red (rows 11 and 12
    Categorical12Step = [1	0.75	0.5 %Light Orange
        1	0.5	0                           % Dark Orange
        1	1	0.6                         % Light Yellow
        1	1	0.2 
        0.7	1	0.55                        % Light Green
        0.2	1	0
        0.65	0.93	1                   % Light Blue
        0.1	0.7	1
        0.8	0.75	1                       % Light Purple
        0.4	0.3	1
        1	0.6	0.75                        % Light Red
        0.9	0.1	0.2];                       % Dark Red

    % Stepped-Sequential
    % Rows 1 to 5 are dark to light red
    % Rows 6 to 10 are dark to light brown
    % Rows 11 to 15 are dark to light green
    % Rows 16 to 20 are dark to light blue
    % Rows 21 to 25 are dark to light purple
    
    SteppedSequential = [ 0.6	0.06	0.06
        0.7	0.175	0.175
        0.8	0.32	0.32
        0.9	0.495	0.495
        1	0.7	0.7
        0.6	0.33	0.06
        0.7	0.438	0.175
        0.8	0.56	0.32
        0.9	0.697	0.495
        1	0.85	0.7
        0.42	0.6	0.06
        0.525	0.7	0.175
        0.64	0.8	0.32
        0.765	0.9	0.495
        0.9	1	0.7
        0.06	0.42	0.6
        0.175	0.525	0.7
        0.32	0.64	0.8
        0.495	0.765	0.9
        0.7	0.9	1
        0.15	0.06	0.6
        0.262	0.175	0.7
        0.4	0.32	0.8
        0.562	0.495	0.9
        0.75	0.7	1];

    %Blue to Grey 8 Step
    %Row 1 is Dark blue
    % Rows 4 and 5 are light blue and near white
    %Row 8 is Dark Grey
    BlueToGrey8Step = [0	0.6	0.8
        0.4	0.9	1
        0.6	1	1
        0.8	1	1
        0.9	0.9	0.9
        0.6	0.6	0.6
        0.4	0.4	0.4
        0.2	0.2	0.2];
    
   %Light  To Dark Blue, 10 Steps
   %Row 1 is Light blue (near white)
   % Row 10 is dark blue 
   LightToDarkBlue10Step = [0.9	1	1
        0.8	0.983	1
        0.7	0.95	1
        0.6	0.9	1
        0.5	0.833	1
        0.4	0.75	1
        0.3	0.65	1
        0.2	0.533	1
        0.1	0.4	1
        0	0.25	1];

    %Blue to Red 18 Step
    %Row 1 is Dark blue
    %Row 18 is Dark Red
    BlueToRed18Step = [0.142	0	0.85
        0.097	0.112	0.97
        0.16	0.342	1
        0.24	0.531	1
        0.34	0.692	1
        0.46	0.829	1
        0.6	0.92	1
        0.74	0.978	1
        0.92	1	1
        1	1	0.92
        1	0.948	0.74
        1	0.84	0.6
        1	0.676	0.46
        1	0.472	0.34
        1	0.24	0.24
        0.97	0.155	0.21
        0.85	0.085	0.187
        0.65	0	0.13];

   %Light To Dark Blue - 7 step
   %Row 1 is Light blue (near white)
   % Row 10 is dark blue 
   LightToDarkBlue7Step = [1	1	1
        0.8	0.993	1
        0.6	0.973	1
        0.4	0.94	1
        0.2	0.893	1
        0	0.667	0.8
        0	0.48	0.6];
    
    %% Assign and return the appropriate matrix(es)
    RampsAll = {GreenToMagenta16Step,Categorical12Step,SteppedSequential,BlueToGrey8Step,LightToDarkBlue10Step,BlueToRed18Step,LightToDarkBlue7Step};
     
    if (nargin==1) && strcmpi('all',lower(varargin{1}))
        varargout = RampsAll{:};
    else
        RampNamesLower = lower(RampNames);
        Indexes = 1:length(RampNames);
        for i=1:nargin
            blIs = strcmp(RampNamesLower,lower(varargin{i}));
            if sum(blIs)==0
                strRampNames = sprintf('%s, ',RampNames{:});
                error('OSUColorRamps:Input Error. Did not find %s. Allowable Inputs are %sor %s', varargin{i},strRampNames,'all');
                varargout = {};
            else
                j = Indexes(blIs);
                varargout{i} = RampsAll{j};
            end
        end
    end
end

