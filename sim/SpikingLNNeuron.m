classdef SpikingLNNeuron < LNNeuron
    properties (SetAccess = private, GetAccess = private)
        inputs;
        kernel;
        nonlinearity;
    end
    
    properties (SetAccess = private)
        membranePotential;
    end
    
    methods
    end
end