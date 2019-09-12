classdef PhaseSaturations < PrimaryVariable
    properties
    end
    
    methods
        function prop = PhaseSaturations(model, varargin)
            prop = prop@PrimaryVariable(model, varargin{:});
            prop.dofname = 's';
        end
    end
end