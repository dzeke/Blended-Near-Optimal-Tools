Interactive-Parallel-Coordinate-Plot
====================================

Makes a parallel coordinate plot and creates a graphical user interface (GUI) to allow a user to interact with the plot.

The major functionalities are:
%
%  1) Option to link objective and decision spaces for an optimization problem and plot all on the same
%       parralel coordinate plot. (Objective Functions and their values [the
%       matrix mObjs] are plotted on the first 1..nO axes while
%       decision variable values in the matrix mDecisions are plotted on nD
%       subsequent axes to the right; pass an empty matrix for mObjs (mObjs=[]) to only show the mDecisions
%       matrix more akin to the standard parallelcoords function).
%
%       A row in each matrix mObjs and mDecisions represents a solution to the underlying optimization
%       problem (i.e., a "point" in the multivariate space).
%
%   2) Tools on the 'Display' tab (far right) to brush, pivot, and interactively manipulate the objective
%       function and decision variable data matricies as well as rows within them (i.e., both by columns/axes and
%       rows/groups/color).
%
%   3) Tools on the 'Interact' tab to generate new solutions or
%       alterantives, e.g., for near-optimal analysis. There are three data sources available to use to generate
%       new solutions: i) the data matrix itself (querying), 2)
%       MATLAB (for linear programs using data specified in the objective
%       function vector(s) and constraint matrix), and 3) executing a GAMS
%       file (note, you will need to install GAMS, see www.gams.com). The
%       methods can generate: a) one solution, b) the maximum extents
%       representing the minimum and maximum values given current selected
%       variables and their values, c) randomly generated solutions, and d)
%       all solutions via enumeration (for MIP problems)
%
%   4) Tools on the 'File' tab to save the current figure settings to the Matlab base workspace and use
%       to recreate the figure
% 
%   5) A variety of parameters/settings to control how the parallel
%       coordinate plot is displayed and labeled (varargin)

% Written in Matlab
%
