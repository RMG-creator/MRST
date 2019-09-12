classdef Pressure < PrimaryVariable
    properties
    end
    
    methods
        function prop = Pressure(model, varargin)
            prop         = prop@PrimaryVariable(model, varargin{:});
            prop.dofname = 'pressure';
        end
    end
end