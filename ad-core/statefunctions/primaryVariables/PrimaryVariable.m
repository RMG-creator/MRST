classdef PrimaryVariable < StateFunction
   
    properties
        dofname
    end
    
    methods
        function prop = PrimaryVariable(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
        end
        
        function value = evaluateOnDomain(prop, model, state)
            % Given state, evaluate the canonical representation for the
            % current model.
            dof   = state.(prop.dofname);
            value = model.disc.evaluateProp(state, dof);
        end
    end
    
end