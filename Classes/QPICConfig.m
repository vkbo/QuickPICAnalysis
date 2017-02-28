
%
%  Class Object :: Holds the QuickPIC config file
% ************************************************
%

classdef QPICConfig
    
    %
    % Properties
    %
    
    properties(GetAccess='public', SetAccess='public')

        Path       = '';     % Path to data directory
        File       = '';     % Config file within data directory
        SimN0      = 1.0e20; % Plasma density for simulation
        PhysN0     = 0.0;    % Plasma density for physics
        Silent     = false;  % Set to true to disable command window output

    end % properties

    properties(GetAccess='public', SetAccess='private')
    
        Name       = '';     % Name of the loaded dataset
        Details    = {};     % Description of simulation. First lines of comments.
        Input      = {};     % Parsed input file
        Constants  = {};     % Constants
        Convert    = {};     % Unit conversion factors
        Simulation = {};     % Simulation variables
        EMFields   = {};     % Electro-magnetic field variables
        Particles  = {};     % Particle variables

    end % properties

    properties(GetAccess='private', SetAccess='private')

        Files      = {};     % Holds possible config files
        Translate  = {};     % Container for Variables class
        NameLists  = {};     % Container for Fortran namelists
        
    end % properties
    
    methods
    end
    
end % classdef

