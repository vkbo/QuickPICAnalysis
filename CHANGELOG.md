# QuickPIC Analysis Toolbox
MATLAB package for analysing QuickPIC data<br>
Current version: 0.1

**Developed by:**<br>
Veronica K. Berglyd Olsen<br>
Department of Physics<br>
University of Oslo

## Development History

### Version 0.1

* Core classes QPICData and QPICConfig added. Derived from the similar classes written for Osiris in the project [OsirisAnalysis](https://github.com/Jadzia626/OsirisAnalysis).
* The QPICData class wraps the dataset and creates a list of what data is available. The data can be extracted by the Data() function for all types of data.
* The QPICConfig class wraps the input file. The class is created by the QPICData class automatically, and its values are contained in the Config struct under the data class. The class assumes per now that the input file is named rpinput, which is the standard name. The class parses all namelists and does an on the fly claculation on conversion constants and makes lists of what diagnostics is available.
