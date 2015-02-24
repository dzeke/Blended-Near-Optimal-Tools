Documentation-readme.txt

The Documentation folder contains the script files and directions to generate each figure in the version of the paper re-submitted in February 2015. Additionally, the folder has PDF versions of the re-submitted and original manuscripts, reviewer comments, and response letter to reviewer comments.

Rosenberg-2015Feb-BlendedNearOptimalTools.pdf : revised version of the manuscript re-submitted to Water Resources Research in Feb 2015 that addresses a 3nd round of reviewer comments.

Rosenberg-2014Dec-BlendedNearOptimalTools.pdf : revised version of the manuscript re-submitted to Water Resources Research in Dec 2014 that addresses 2nd round of reviewer comments.

Rosenberg-2014Aug-BlendedNearOptimalTools.pdf : revised version of the manuscript re-submitted to Water Resources Research in August 2014 that addresses 1st round of reviewer comments.

Rosenberg-2013Aug-BlendedNearOptimalTools.pdf : original version of the manuscript submitted in August 2013 to Water Resources Research.

Rosenberg-2015Feb-ResponseLetter3-BlendedNearOptimalTools.pdf : Letter listing reviewer comments, author responses and descriptions of changes made in the Feb, 2015 version of the paper to address the comments.

Rosenberg-2014Dec-ResponseLetter2-BlendedNearOptimalTools.pdf : Letter listing reviewer comments, author responses and descriptions of changes made in the Dec, 2014 version of the paper to address the comments.

Rosenberg-2014Aug-ResponseLetter1-BlendedNearOptimalTools.pdf : Letter listing reviewer comments, author responses and descriptions of changes made in the August, 2014 version of the paper to address the comments.


LoadYourOwnModel.m : Directions and example for how to load your own model data into the near-optimal tools. Specific directions for linear programs, more general directions for mixed-integer and non-linear programs.


Within the ScriptsForPaper folder:

- Fig_GenForNearOptPaper.m : Matlab script to use to generate Figures 1, 2, 3, and 5 in the revised paper. To run this script, you must download all files in the AlternativeGeneration, InteractiveParallelPlot, and EchoReservoirApplication folders. In Matlab, add the folder where you downloaded the files to your Matlab path, set the Matlab directory to the same location, and enter the following command at the Matlab command prompt:

	>> FigGenForNearOptPaper

The file also has additional directions to interactively generate Figures 4, 5, 6, and 7 from Figure 3 as comments in the .m file. Figure 5 can also be automatically generated.

The file also contains the commands and directions to generate all the results that are discussed in the manuscript but not presented.

- Fig1_FeasibleNearOptCompare.m : Matlab script to generate Figures 1 and 2 in the revised paper. See Fig_GenForNearOptPaper.m for the parameter settings used for the paper figure.

- doMGA.m : a Matlab file that includes the logic for generating alternatives by various MGA methods. Used to generate Figures 1, 3, 4, 5, and 6.

- Delcols.m, extrdir.m, extrpts.m, polygeom.m : other Matlab files used by the above scripts.

