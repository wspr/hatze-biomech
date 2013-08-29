# A Matlab implementation of Hatze's 1980 anthropometric body segment parameter model

Biomechanics simulations often require accurate models of parameters such as mass, moment of inertia, etc., for individual body segments of an individual.
H. Hatze published an anthropometric model of the human body in 1979/1980 to estimate these body segment parameters.
His model was composed of 17 body segments with varying cross-section and density.
It accurately models the masses, centroids, moments of inertia, products of inertia (where applicable), and principle axes (where applicable) for these body segments using a total of 242 anthropometric measurements of an individual.
It is the most detailed such model developed, but is little used in the research community due to the large number of measurements required and the mathematical complexity involved in using those measurements.

The code herein replicates (or will) Hatze's model in Matlab, providing a public implementation of his model for the first time.
We hope that this will increase the viability of using his model.

# Code layout

The code in the `master` branch is currently at a state which can draw all of the geometry of a set of measurements (see `person_plot.m`) and intertwined are many of the calculations required to obtain volume, mass, moment of inertia, etc.
These calculations can be interpolated in the cases of the limbs (see MODSIM 2013 paper) to achieve better accuracy especially in the cases with fewer measurements.
Many other calculations, however are still incomplete or, worse, incorrect; this shall hopefully be rectified in the coming weeks (time of writing Sept 2013).

Since submission of the MODSIM paper, changes have been made to allow joint angles to be reflected in the visual output.
This functionality is still half-baked, however; tweaking the drawing functions to display correctly at arbitrary orientations is still ongoing.

## Manifest

* `hatze-cheatsheet.pdf` — a short summary taken from Hatze's works
* `hatze_meas.txt` — transcription of Hatze's measurements
* `person_generate.m` — main function
* `person_*.pdf` — driver files to:
   * `person_plot.m` — generic ‘run this file to see the results’
   * `person_hatze_compare.m` — compare our calculations to Hatze's
   * `person_reduce.m` — generate the results for the MODSIM paper
   * `person_move.m` — experiment with joint angles
* `segment_*.m` — the functions used to calculate and plot the various body segments
* `plot_*.m` — various plotting functions for the graphical output

Miscellaneous functions

* `rotation_matrix_zyx.m` — XYZ rotation matrix

# TODO

* Complete all segment parameter calculations (ongoing)
* Allow joint angle inputs (ongoing)
* Develop a generalised extruded ellipse shape to avoid interpolation and improve accuracy
* Verify with measured data
* Incorporate centre of mass coordinate systems (HOMSIM report, Hatze 1981)
* Provide reference documentation for all equations and shape derivations
* Compare against models of Hanavan, Yeadon, etc.
* Provide a method for taking measurements from photographic (and/or rangefinder?) data
* Consider adding additional degrees of freedom (e.g., knee ab/adduction)
* Add knees
* Improve feet (Dillon 2001)
* Add separate fingers (e.g., see new Optitrack software with digit motion tracking)
* Improve head (notably jaw)

And of course:

* Turn it into a dynamic model

# Authors and copyright

This work is freely distributed under the terms and conditions of the Apache Licence, v2:
 <http://www.apache.org/licenses/LICENSE-2.0>

Jianan (Tony) Wang transcribed many of Hatze's equations for this work and wrote initial code to do some of the plotting.
Will Robertson wrote most of the code relating to the graphics objects and structured the program into its current form.

# References

* H. Hatze (1979) ‘A model for the computational determination of parameter values of anthropomorphic segments’, CSIR Techn. Report TWISK 79, Pretoria

* H. Hatze (1980) ‘A mathematical model for the computational determination of parameter values of anthropomorphic segments’, Journal of Biomechanics, vol. 13, pp. 833-843

* H. Hatze (1981) ‘HOMSIM: a simulator of three-dimensional hominoid dynamics’, CSIR Techn. Report SWISK 23, Pretoria

* M. Dillon (2001) ‘Biomechanical Models for the Analysis of Partial Foot Amputee Gait’ PhD Thesis, Queensland University of Technology, Australia
