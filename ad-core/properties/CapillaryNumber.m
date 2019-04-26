classdef CapillaryNumber < GridProperty
    
    properties
    end
    
    methods
        function gp = CapillaryNumber(varargin)
            gp@GridProperty(varargin{:});
        end
        function Nc = evaluateOnDomain(prop, model, state)
        % For this grid property, it is expected that it has been set before.
            Nc = state.CapillaryNumber;
        end
    end
end
