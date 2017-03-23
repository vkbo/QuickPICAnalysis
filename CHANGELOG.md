# QuickPIC Analysis Toolbox
MATLAB package for analysing QuickPIC data<br>
Current version: Dev0.3

**Developed by:**<br>
Veronica K. Berglyd Olsen<br>
Department of Physics<br>
University of Oslo

## Development History

### Version 0.3
_Under Development_

### Version 0.2
_Released 23.03.2017_

* Added QPICTools class that is essentially a namespace containing static functions necessary for various QuickPIC related operations like labelling axes, etc.
* Added QPICType superclass that handles conversions, slicing, lineouts and other common functions on QuickPIC datasets. This class is extended by three subclasses:
    * QPICScalar: A class that extends QPICType for scalar type data like charge and potential.
    * QPICVector: A class that extends QPICType for vector type data like fields and current.
    * QPICBeam: A class that extends QPICType for raw beam dump.
* Added standard plots for vector and scalar data.

### Version 0.1
_Released 07.03.2017_

* Core classes QPICData and QPICConfig added. Derived from the similar classes written for Osiris in the project [OsirisAnalysis](https://github.com/Jadzia626/OsirisAnalysis).
* The QPICData class wraps the dataset and creates a list of what data is available. The data can be extracted by the Data() function for all types of data.
* The QPICConfig class wraps the input file. The class is created by the QPICData class automatically, and its values are contained in the Config struct under the data class. The class assumes per now that the input file is named rpinput, which is the standard name. The class parses all namelists and does an on the fly calculation on conversion constants and makes lists of what diagnostics is available.
