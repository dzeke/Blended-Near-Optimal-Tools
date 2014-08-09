GenerateAlternatives-readme.txt

The GenerateAlternatives folder contains the Matlab scripts needed to generate alternatives by the stratified, Monte-Carlo Markov Chain (MCMC) method for a closed, bounded region defined by a system of linear inequalities using the MCMC Gibbs method. To use, call the function

	 [X,vValid] = stratgibbs(p,A,b,options)

Where p is the desired number of samples, A is the m x n matrix defining coefficients for m constraints on n variables, and b is and m x 1 vector of the right-hand side values for the constraint. Constraints are thus in the form Ax <= b. options is a structure of optional parameters that determine how to sample (see the stratgibbs.m file for further info).

The function returns X as a p x n matrix of the sampled alternatives where each row is an alternative. vValid is a p x 1 vector of codes describing the validity of each sampled alternative including whether the alternative is feasible according to the constraints Ax<=b. A description of the codes for vValid is provided in stratgibbs.m.

stratgibb.m handles the initial stratification across all problem dimensions and dimension reduction. See further documentation on use in stratgibbs.m and on the method in section 3 of the paper. stratgibbs calls:

maxextentind.m to calculate the separate and independent minimum and maximum extents for each variable within the region. See Section 3 step b of the paper.

maxextentgibbs.m to Gibbs sample from the reduced region specified by the stratified values defined in stratgibbs.m.

