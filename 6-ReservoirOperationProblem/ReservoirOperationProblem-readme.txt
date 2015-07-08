ReservoirOperationProblem-readme.txt

The ReservoirOperationProblem folder contains the data, model, and directions to undertake optimal and near-optimal analysis for a simple multipurpose reservoir operation problem involving generating hydropower and irrigation benefits over 6 time periods.

- ReservoirOperationProblem.m : A matlab script and function that defines the system of linear constraints defining the reservoir operation problem, solves for the optimal solution, and plots othe optimal solution in the parallel coordinate plotting tool. Directions in the file 1-Documentation/NearOptimalManagement-Lab.pdf describe how to interactively use the parallel coordinate plotting tool to undertake a variety of near-optimal analyses. At the Matlab prompt, enter the command:

	>> ReservoirOperationProblem([])

to generate the optimal solution and plot from the default data listed in the .pdf file. To run with custom data, substitute a Matlab structure for the [] with fields listed in ReservoirOperationProblem.m (see header of file).
