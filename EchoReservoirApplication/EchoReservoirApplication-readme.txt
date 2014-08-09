EchoReservoirApplication-readme.txt

This folder contains the data, model, and script files used for the water quality management application for Echo Reservoir that illustrates use of the blended near-optimal alternative generation, visualization, and interaction tools. These files are:

WaQualityModelGLOB_NE5.gms – a General Algebraic Modeling System (GAMS) file that has all the input data, model equations that define the original single-objective model used to find the optimal solution, model to produce the Modeling to Generate Alternatives, model for the multi-objective formulation, and loops used to generate the MGA and pareto alternatives. The GAMS file output is the file WQNE_outG6.gdx.

- WQNE_outG6.gdx : GAMS file with all the model input data and model results for the optimal solution, MGA and pareto alternatives.

- LoadEchoGamsResults.m : Matlab script  that reads the GAMS results in WQNE_outG6.gdx, generates additional alternatives using the stratified sampling tool (see Alternative Generation), and compiles and formats all of the results for plotting in the interactive parallel plotting tool. Part of the formatting is turning all decision variable values into units of phosphorus removed for more convenient plotting on the parallel coordinate plot. To load the results and generate Figure 3 in the paper, make sure the folder where you downloaded WQNE_outG6.gdx is the active Matlab folder. Then enter:

       >> LoadEchoGamsResults('WQNE_outG6.gdx',3,2500,0) 

Here 3 is the number of the optimal solution in the .gdx file, 2500 alternatives are sampled, and 0 is a flag to read and format for the single-objective formulation that only minimizes reduction cost.

To load the results and generate Figure 5 in the paper for the multi-objective problem, enter:

       >> LoadEchoGamsResults('WQNE_outG6.gdx',3,2500,1)

To use these tools, download all files in the AlternativeGeneration and InteractiveParallelPlot folders.

