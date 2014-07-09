Interactive-Parallel-Coordinate-Plot
====================================
Makes a parallel coordinate plot and creates a graphical user interface (GUI) to allow a user to interact with the plot to rearrange axes, highlight rows, and generate new alternatives.

The major functionalities are:

1) Option to link objective and decision spaces for solutions to an optimization problem and plot all on the same parralel coordinate plot. (nO number of Objective Functions and their values are plotted on the first 1..nO axes while nD number of decision variables and their values are plotted on the subsequent axes to the right).

2) Tools to control the displaying, brush, pivot, and interactively manipulate the objective function and decision variable data matricies as well as rows within them (i.e., both by columns/axes and rows/groups/color). 

3) Tools to select particular solutions (or attributes of solutions) and generate new alterantives, e.g., for near-optimal analysis.

4) Tools to save the current figure settings to the Matlab base workspace and use to recreate the figure

5) A variety of parameters/settings to control how the parallel coordinate plot is displayed and labeled (varargin)

Written in Matlab

