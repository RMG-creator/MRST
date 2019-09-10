classdef Pressure < StateFunction
   
    properties
        name = 'pressure';
    end
    
    methods
        function prop = Pressure(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
        end
    end
end