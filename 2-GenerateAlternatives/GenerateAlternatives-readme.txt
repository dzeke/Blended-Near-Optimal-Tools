GenerateAlternatives-readme.txt

The GenerateAlternatives folder contains the Matlab scripts needed to generate alternatives by the stratified, Monte-Carlo Markov Chain (MCMC) method for a closed, bounded region defined by a system of linear inequalities using MCMC sampler. To use, download all the files in the folder to a folder on your local machine. Add the target folder to your Matlab path, then, at the Matlab prompt, enter:

	>> [X,vValid] = stratgibbs(p,problemstruct,options)

Where p is the desired number of samples, problemstruct is a structure with fields for the inequality, equality, lower-bound, and upper-bound constraints that define the polytope (see documentation in stratgibbs). options is a structure of optional parameters that determine how to sample (see the stratgibbs.m file for further info).

The function returns X as a p x n matrix of the sampled alternatives where each row is an alternative. vValid is a p x 1 vector of codes describing the validity of each sampled alternative including whether the alternative is feasible according to the constraints Ax<=b. A description of the codes for vValid is provided in stratgibbs.m.

stratgibbs.m handles the initial stratification across all problem dimensions and dimension reduction. See further documentation on use in stratgibbs.m and on the method in section 3 of the paper. stratgibbs calls:

maxextentind.m to calculate the separate and independent minimum and maximum extents for each variable within the region. See Section 3 step b of the paper.

A function (file) based on the Monte Carlo Markov Chain sampling method set in options.MCMCMethod (see stratgibbs.m for details). The files are:

   - cprnd.m to Gibbs or Hit-and-run sample

   - maxextentgibbs.m to Gibbs sample by the method of maximum extents when cycling through coordinate directions from the reduced region specified by the stratified values defined in stratgibbs.m.

Additional files to generate alterantives are:

   - chebycenterFull.m : used to find an initial point inside the polytope to start Monte Carlo Markov chain sampling.

   - OptimiFull.m : used to check that the fields in problem structure (problemstruct) are valid.
