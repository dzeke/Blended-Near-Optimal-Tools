Blended Near-Optimal Tools
==============================
This repository stores the Matlab 2013a source code for blended near-optimal tools that (i) generate alternatives, (ii) visualize alternatives and allow a user to interactively explore the region from which the alternatives are generated. The repository also contains the (iii) data and model files for an example water quality management application for Echo Reservoir, Utah (linear program) as well as further (iv) documentation (paper submitted to Water Resources Research, scripts and directions to generate each figure in the paper that describes these tools).

1) GENERATE ALTERNATIVES

Uses stratified Monte-Carlo Markov Chain sampling to identify a large number of near-optimal alternatives that comprehensively span the region through both the decision and objective spaces.

2) PARALLEL COORDINATE VISUALIZATION and INTERACTION

Parallel coordinate plot places axes for all objectives and decision variables side-by-side on a single page and shows the generated alternatives across the decision and objective spaces. Interaction controls on and next to the plot:
-	Render generated alternatives (mouse-over to read a value, highlight individual or groups of alternatives on the plot, and specify the axes order from left to right across the plot).
-	Allow the users to explore the near-optimal region (set sliders on axes to specify features, dynamically update the model formulation, and generate individual or families of new alternatives with specified features) 
-	Direct exploration into different parts of the region (re-order axes or relax the constraint that specifies the size of the region.
-	Save figure settings to the Matlab e workspace and use to recreate the figure
-	Optional parameters/settings to control how the plot is displayed and labeled (varargin)
Tools are described here for near-optimal optimization problems but are general and can be used for any high-dimensional data set or system of equations.

3) ECHO RESERVOIR WATER QUALITY MANAGEMENT APPLICATION

Illustrates use of the tools for a linear program to identify the cost-effective phosphorus removal practices to reduce the phosphorus load to Echo Reservoir in the Weber basin, Utah to a level specified in a pending Total Maximum Daily Load (TMDL) program for the reservoir. Includes all the data and model files and the script that moves data from the models to the near-optimal tools.

4) DOCUMMENTATION (PAPER and SCRIPTS)

Includes the paper submitted to Water Resources Research (under review) that describes this work as well as peer-reviewer comments and author responses. Also provides scripts and directions to generate each figure shown in the paper. Further documentation on each tool is provided in the sub-directory of each tool (XXXX-readme.txt) as well as at the top of each source code file.

==================
CITATION

David E. Rosenberg (in review) "Near-optimal alternative generation, visualization, and interaction for water resources decision making". Water Resources Research. Submitted August 2014.

LICENSING

All code is distributed AS-IS with no expressed or implied warranty regarding functionality. The code or parts may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires prior written permission from the author.

BUG REPORTS and FEEDBACK

This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker for this repository. And note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.

David E. Rosenberg
Department of Civil & Environmental Engineering and Utah Water Research Lab
Utah State University
Email: david.rosenberg@usu.edu


