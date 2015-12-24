InteractiveParallelPlot-readme.txt

The InteractiveParallelPlot folder contains the Matlab functions to create a Parallel Coordinate plot that places axes for all objectives and decision variables for an optimization problem side-by-side on a single page. Additionally, the function places interactive controls on and next to the plot to allow the user to control the plot rendering, explore the region, update the model formulation, and generate new alternatives from the updated formulation. To use, download all the files in this folder as well as the GenerateAlternatives folder. Than, at the Matlab prompt, call the main function nearoptplotmo2.m as follows:

	 >> [hWindow] = nearoptplotmo2(mObjs,mDecions,OptionList)

where mObjs is a p by nO matrix of objective function values for the nO objectives and p alternatives (a row is an alternative), mDecisions is a p by n matrix of n decision variable values for the p alternatives (together, a row in mObjs and in mDecisions define the objective and decision variable components of an alternative). OptionList is an optional cell array of paired items {‘parameter1’,value1,’parameter2’,value2, … } that control numerous plotting, rending, interaction, and other options. See full documentation in the header of file nearoptplotmo2.m. These options include providing data that define the problem and near optimal region (constraint set, objective functions, tolerance, etc…) which the figure can then use to streamline the process to update the model formulation and generate new alternatives.
To use this tool, you will also need to download all the files in the GenerateAlternatives folder.

Interaction controls are organized into six categories on the plot:

A) Plot Data menu. Click Save Data (Cntrl-D) to save all data and current plot settings to the base workspace. The command window will also show the command to enter to re-create the figure. Click Load Data... (Cntrl-L) to load data for additional alternatives from the base workspace into the figure. Click Print to PDF... to print a high-resolution version of the current figure to a PDF file you specify.

B) Controls Menu. Menu items to show/hide various plot features and controls including the control panel at right, sliders, checkboxes, and grouping table at bottom.

C) Update Formulation Menu. Menu items to streamline the process to add one or more new optimization model objective functions or constraints and generate new alternatives from the updated formulation (must provide description of existing model formulation within the OptionList, see above).

D) Display tab. Controls the  rendering and formatting of plotted data including

	o Font size – font size for text labels

	o Move/Delete Axes. Check one or more checkboxes below axes on the plot and then use the controls to move, reorder (either by categories specified in optional data provided in mActCat [see documentation in nearoptplotmo2.m] or by the ranges of values for each axis) or delete axes.

	o Group Traces to control the order and rending of groups of traces. Groups are specified by the optional parameter vGroup (a p x 1 vector). The group with the highest Order # is plotted last (on top). Use the Highlight Group selection box to assign a group to highlight (plot last in thick black lines).

E) Color Ramp tab. Add a color ramp to one axis to facilitate comparing values between that axis and each other axis.

	o On the plot, check the axis for which to color ramp values.

	o Then on the Color Ramp tab, set the number of color bands to include in the ramp

	o Specify the direction of the ramp (Ascend = ascending values get progressively darker, e.g., for an objective to maximize; Descend = descending values get progressively darker, e.g., for an objective to minimize). Darker colors are plotted last on top andi indicate closer to optimal.

	o Check the Ramp color box. The traces will re-color, a color ramp will appear on the checked axis (colors on other axes will vary depending on the relationship with the selected axis). Below the color ramp controls, new checkboxes will appear that define the ranges for each color band in the ramp (both as an absolute value and % of the value associated with the darkest band). Uncheck a box to hide traces for the color band.

F) Interact tab. Controls to explore the near-optimal region shown on the plot and further defined by the optional data:
    - For linear programs in Matlab: AMat, Brhs, cFunc, NearOptConstraint, and OptSolRow provided as input to the function in OptionList. The first three parameters define the model constraints and objective functions while the 4th and 5th paramters identify the near-optimal tolerance constraint within the constraint set and row in the data set that is the optimal solution (see documentation in nearoptplotmo2.m).
    - For models written in the General Algebraic Modeling System (GAMS) : sGamsFile (name of the gams file with the model)

The interactive tools are:

	o Interactive exploration with sliders. Turn on the sliders (Controls menu=>Sliders=>uncheck Hide sliders). The sliders will appear and each height shows the allowable range for the variable or objective. Set the value for one slider within it's allowable range. The checkbox for the axis will check to indicate the slider is "set". Other slider heights will update. Set the slider for a second variable within it's range and continue until setting some or all of the sliders. Uncheck a box below an axis to unset the slider (ranges for sliders will subsequently adjust). Alternatively, in the Filter Existing Alternatives box at the right, enter a "Set value" for an axis and the slider will adjust to that value.

	o Near Optimal Tolerance. Enter a new tolerance value (gamma in the paper), then in the Generate New Alternatives box, click the Generate button to sample new alternatives from the region defined by the entered tolerance value. Note, for multi-objective problems, the near-optimal tolerance for each objective is instead entered in text boxes that appear in the Filter Existing Alternatives box to the right of the listing of objectives. 

	o Subspace Error. An error applied to new constraints added to the model formulation by setting slider values (see below). A value of zero adds an equity constraint (decision variable value = slider setting) whereas a Subspace Error > 0 adds two constraints that limit decision variable values to the slider setting +/- (Subspace Error)*(Tick Spacing) [error is meausred as a fraction of the tick spacing).

	o Generate New Alternatives. Controls allow the user to specify:
		- Type of generation as either a single alternative, maximum extents, random sample, or enumerate all (only for mixed-integer problems)
		- Engine to use to generate alternatives as either from existing data on the plot, Matlab (for linear problems) using Matlab's linprog, or the General Algebraic Modeling System (GAMS).
		- For random sampling, specify the # samples (number of alternatives to generate)
		- If using GAMS, specify the full path of the GAMS file (note, the file will need to be organized in a particular way so that Matlab can pass on input data. See further documentation in the Matlab file EnumNEIntSolsGams4.m and the example in the AmmanJordanApplication folder.
		- Generate button. Will generate the new samples according to the slider settings, checked axes, generation type and engine.

	o Filter Existing Alternatives. Show either:

		- Individual alternative (row in the input data). On the Controls menu=>select Alternatives=>Show current. Then use the arrows to scroll through the alternatives. Or enter the row number in the Set to box. The trace will highlight on the plot and values show in the Set value text boxes.

		- Show filtered (all alternatives meeting criteria the criteria specified by the slider settings, checked axes, and SubspaceError). Turn on from the Controls menu=>select Alternatives=>Show filtered. Groups of alternatives already on plot that meet the criteria will highlight in a dark color. Click the generate button in the Generate New Alternatives box to generate additional alternatives meeting these same criteria.


Finally, mouse over any axis to read the value in real units.

Again, all of the above interactive features can also be supplied as optional parameters in OptionList so that the plot opens with the feature active. See further documentation of OptionList within the header of the file nearoptplotmo2.m.

ExampleTestPlots.m -- contains multiple commands that generate parallel coordinate plots that use several of the above features.

All the other files in this folder are called by nearoptplotmo2.m to produce a plot. See further descriptions in the file headers.
