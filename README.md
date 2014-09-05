Blended Near-Optimal Tools
==============================
This repository stores the Matlab 2013a source code for blended near-optimal tools that (a) generate alternatives, (b) visualize alternatives, and allow a user to interactively explore the near-optimal region from which alternatives are generated. The repository also contains the data and model files for a (c) linear programming example application to manage water quality for Echo Reservoir, Utah and (d) mixed-integer programming example application to manage water supply and demands in Amman, Jordan. The repository additionally provides (e) documentation on the above tools--including the paper submitted to Water Resources Research, scripts, and directions to generate each figure in the paper.

Near-optimal alternatives perform within a (near-optimal) tolerable deviation of the optimal objective function value and are of interest to managers and decision makers because near-optimal alternatives can address un-modelled objectives, preferences, limits, uncertainties, or issues that are not considered by the original optimization model or it's optimal solution. Mathematically, the region of near-optimal alternatives is defined by the constraints of the original optimization model as well as a constraint that limits alternatives to those with objective function values that are within a tolerable deviation of the optimal objective function value. The code and tools within this repository allow users to generate and visualize the structure and full extent of the near-optimal region to an optimization problem. The tools also allow users to interactively explore region features of most interest, streamline the process to elicit un-modelled issues, and update the model formulation with new information. The tools and their use are described here for generating, visualizing, and interactively exploring near-optimal alternatives to optimization problems, but the tools are general and can be used to generate and visualize points within any high-dimensional, closed, bounded region that can be defined by a system of constraints. The parallel coordinate visualization and several interaction tools can also be used for any high-dimensional data set.   

==============================
1) DOCUMMENTATION (Paper and Scripts)

Includes the paper submitted to Water Resources Research (under review) that describes this work as well as peer-reviewer comments and author responses. Also provides scripts and directions to generate each figure shown in the paper. These scripts make use of each of the tools described in sections 2) to 4). To generate each figure in the manuscript:

    a. Download all the files in the folder 1-Documentation\ScriptsForPaper as well as the folders 2-GenerateAlternatives, 3-InteractiveParallelPlot, and 4-EchoReservoirApplication folders to a single folder on your local machine.
    b. In Matlab, add the target folder where you downloaded the files to your Matlab path
    c. Set the Matlab directory to the same folder
    d. At the Matlab command prompt, enter:

	       >> FigGenForNearOptPaper

    e. The script will generate Figures 1, 2, 3, and 6 in the submitted manuscript.
    f. See the file Documentation-readme.txt for further instructions on how to interactively generate figures 4, 5, and 7 from Figure 3 and 6.

==============================
2) GENERATE ALTERNATIVES

Uses stratified Monte-Carlo Markov Chain sampling to identify a large number of alternatives that comprehensively span the near-optimal region through both the decision and objective spaces. To get started:

    a. Download all the files in the folder 2-GenerateAlternatives to a folder on your local machine.
    b. Add the target folder to your Matlab path.
    c. At the Matlab prompt, invoke the function
       
         >> stratgibbs(...)
       
    d. See additional directions in the file GenerateAlternatives-readme.txt for use of the stratgibbs function.

==============================
3) PARALLEL COORDINATE VISUALIZATION and INTERACTION

Parallel coordinate plot places axes for all objectives and decision variables side-by-side on a single page and shows the generated alternatives across the decision and objective spaces. Interaction controls on and next to the plot:
-	Render generated alternatives (mouse-over to read a value, highlight individual or groups of alternatives on the plot, and specify the axes order from left to right across the plot).
-	Allow the user to explore the near-optimal region (set sliders on axes to specify features, dynamically update the model formulation, and generate individual or families of new alternatives with specified features) 
-	Direct exploration into different parts of the region (re-order axes or relax the constraint that specifies the size of the region.
-	Save figure settings to the Matlab e workspace and use to recreate the figure
-	Optional parameters/settings to control how the plot is displayed and labeled (varargin)

To get started:

    a. Download all the files in the folders 2-GenerateAlternatives and 3-InteractiveParallelPlot to a single folder on your local machine.
    b. In Matlab, add the target folder where you downloaded the files to your Matlab path
    c. At the Matlab command prompt, invoke the following function:

	       >> nearoptplotmo2(...)

    d. See additional directions in the file InteractiveParallelPlot-readme.txt and in the header of the nearoptplotmo2.m for use of the function.

==============================
4) ECHO RESERVOIR APPLICATION (Linear program to manage water quality)

Illustrates use of the tools for a linear program to identify the cost-effective phosphorus removal practices to reduce the phosphorus load to Echo Reservoir in the Weber basin, Utah to a level specified in a pending Total Maximum Daily Load (TMDL) program for the reservoir. Includes all the data and model files and the script that moves data from the models to the near-optimal tools.

To use:

    a. Download all the files in the folders 2-GenerateAlternatives, 3-InteractiveParallelPlot, and 4-EchoReservoirApplication to a single folder on your local machine.
    b. In Matlab, add the target folder where you downloaded the files to your Matlab path.
    c. Set the Matlab directory to the same folder
    d. At the Matlab command prompt, enter:

	       >> LoadEchoGamsResults('WQNE_outG6.gdx',3,2500,0)

    e. Matlab will read optimal and modelling to generate alternatives (MGA) results from the file WQNE_outG6.gdx. Matlab will also generate 2,500 near-optimal alternatives (see Section 2 above) and plot the optimal, MGA, and near-optimal results in the Interactive Parallel Plotting tool (see Section 3 above).
    f. See additional directions in the file EchoReservoirApplication-readme.txt.

==============================    
5) AMMAN, JORDAN APPLICATION (Mixed-integer program to manage water supply/demand)

Illustrates use of the tools for a mixed-integer program to identify the cost-effective combination of new water supply and conservation strategies to balance water supplies and demands in Amman, Jordan through 2020. Includes all the data and model files and scripts to move data from the models to the near-optimal tools and to run the model from within the Interactive Parallel Coordinate Plotting tool.

To use:

    a. Download all the files in the folders 3-InteractiveParallelPlot and 5-AmmanJordanApplication to a single folder on your local machine.
    b. In Matlab, add the target folder where you downloaded the files to your Matlab path.
    c. Set the Matlab directory to the same folder
    d. At the Matlab command prompt, enter:

	       >> LoadAmmanJordan(AmmanJordanUtilOptNear.gdx)

    e. Matlab will read optimal and near-optimal results from the file AmmanJordanUtilOptNear.gdx and plot the results in the Interactive Parallel Plotting tool (see Section 3 above).
    f. For interactive use, see additional directions in the file AmmanJordanApplication-readme.txt.

==================
CITATION

David E. Rosenberg (in review) "Near-optimal alternative generation, visualization, and interaction for water resources decision making". Water Resources Research. Submitted August 2014.

LICENSING

All code is distributed AS-IS with no expressed or implied warranty regarding functionality. The code or parts may be used for non-commercial purposes so long as the use is cited per the citation above. Use for any commercial purpose requires prior written permission from the author.

BUG REPORTS and FEEDBACK

This code is possibly laden with bugs so bug reports and feedback are much appreciated. Please submit via the the issue tracker for this repository. And note, that while much appreciated, there is no promise of when--or if--the bug will be fixed.

===================
AUTHOR and CONTACT

    David E. Rosenberg
    Department of Civil & Environmental Engineering and Utah Water Research Lab
    Utah State University
    Email: david.rosenberg@usu.edu


