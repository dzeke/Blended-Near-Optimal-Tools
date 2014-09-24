Documentation-readme.txt

The Documentation folder contains the script files and directions to generate each figure in the version of the paper re-submitted in August 2014. Additionally, the folder has PDF versions of the re-submitted and original manuscripts, reviewer comments, and response letter to reviewer comments.

Rosenberg-NearOptimalWaterResourcesDecisionMaking-Aug2013.pdf : original version of the manuscript submitted in August 2013 to Water Resources Research.

Rosenberg-NearOptimalWaterResourcesDecisionMaking-Aug2014.pdf : revised version of the manuscript re-submitted to Water Resources Research in August 2014 that addresses reviewer comments.

Rosenberg-ResponseLetter-NearOptimal2-Aug2014.pdf : Letter listing reviewer comments, author responses and descritions of changes made in the August, 2014 version of the paper to address the comments.


Within the ScriptsForPaper folder:

- Fig_GenForNearOptPaper.m : Matlab script to use to generate Figures 1, 2, 3, and 5 in the revised paper. To run this script, you must download all files in the AlternativeGeneration, InteractiveParallelPlot, and EchoReservoirApplication folders. In Matlab, add the folder where you downloaded the files to your Matlab path, set the Matlab directory to the same location, and enter the following command at the Matlab command prompt:

	>> FigGenForNearOptPaper

Additional directions to interactively generate Figures 4, 5, 6, and 7 from Figure 3 are listed in comments in the .m file. Figure 5 can also be automatically generated.

- Fig1_FeasibleNearOptCompare.m : Matlab script to generate Figure 1 in the revised paper. See Fig_GenForNearOptPaper.m for the parameter settings used for the paper figure.

- Fig2_CarParCompare.m : Matlab script to generate Figure 2 in the revised paper. See Fig_GenForNearOptPaper.m for the parameter settings used for the paper figure.

- Delcols.m, extrdir.m, extrpts.m, polygeom.m : other Matlab files used by the above scripts.

