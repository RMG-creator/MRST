classdef PhaseSaturation < StateFunction
   
    properties
        name = 's';
    end
    
    methods
        function prop = PhaseSaturation(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
        end
    end
end