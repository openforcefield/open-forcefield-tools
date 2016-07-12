Contains

* thermopyl_description.txt: a description of the data representation of single molecule simulation data for our simple test case

* bonded_term_list.py: take multiple parmtops and generate a pandas dataframe as output by thermopyl

* trajectory_parser.py: take multiple .nc files with identifier names and a pandas dataframe with property names for single atom bonded properties (including the atom numbers) and populate those property pandas dataframe.

* generate_single_molecule_properties.py: sample script running bonded_term_list.py and trajectory_parser.py on a specific data set.

* * Bond equilibrium lengths
* * Bond equilibrium standard deviations
* * Angle equilibrium angles
* * Angle equilibrium standard deviations
* * Torsion Fourier representations (0, 1, 2, 3, 6 + phase angle) (determined by least square fit)

and populate the pandas dataframe.
