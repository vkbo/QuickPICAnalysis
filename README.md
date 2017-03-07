# QuickPICAnalysis

MATLAB code for analysing QuickPIC data.

The Open Source version of QuickPIC is available [here](https://github.com/UCLA-Plasma-Simulation-Group/QuickPIC-OpenSource).

### Developed By

Veronica K. Berglyd Olsen<br>
Department of Physics<br>
University of Oslo

Derived from code from the [OsirisAnalysis](https://github.com/Jadzia626/OsirisAnalysis) project.

## Usage

As OsirisAnalysis, QPICAnalysis is layered. The core classes simply wrap the dataset and inputfile. The second layer contains data processing classes which does unit conversion and statistics, and the third layer is a set of standard plot functions that build on the previous two.

Currently only the first layer is ported from OsirisAnalysis.

### Layer 1

The basic layer of QPICAnalysis. Contains two classes:

#### QPICData

A wrapper class for the simulation data folder. It scans for data dumps (HDF5 files) and builds a tree of the available datasets. The location of the dataset is set with the Path method in the class. This method takes either a full path to the datafolder, or if only the name of the set is specified, will look it up in one of the folders specified in QPICSettings.m.

All data can be extracted from the method Data, which takes as input:

1. Time: The dump number of the data requested.
2. Type: The data type requested. Valid inputs are:
    * F: Electromagnetic fields
    * J: Current density
    * PSI: Potential
    * Q: Charge
    * RAW: Particle dumps
3. Set: Applies to type F and J and denotes which vector component to extract. Valid inputs are B, E or J with X, Y or Z as a second character for vector. E.g. BX, EY, JZ, etc.
4. Species: Specify which species to load. Valid inputs are:
    * EB: For beam density data
    * EPxx: For plasma density where xx denotes plasma number in two digits
    * EBxx: For raw particles where xx denotes which beam number in two digits
5. Slice: If 2D slices are specified, these can be loded for all grid based dumps with either XY, XZ or YZ. For 3D data leave blank.

All fields that do not apply should be specified as empty string, i.e. ''.

For an overview of available datasets, there is a struct named SimData.Data in the class that contains a list of all available sets with path and number of files in each.

#### QPICConfig

This class wraps the rpinput file that should be located in the same folder as the dataset.

An instance of this class is automatically created by the QPICData class when a dataset is loaded. It is available under the class object as Config.

The object parses the namelists contained in the input file and extracts all simulation variables. The majority of these are then parsed and organised into a Simulation, Beam, Plasma, and Diag struct containing the relevant inputs in SI units. These structs will fall back to default values if the setting is not available in the input file, so a lookup should never return a not found error.

A number of conversion constants are also calculated on the fly and are available in the struct named Convert. The key natural constants are also available in the struct Constants.

The raw values from the input file, with its original field name, are contained in the struct Input. This tree is built directly from the input file as is. I.e. there is no guarantee the field exists or holds a default value.
