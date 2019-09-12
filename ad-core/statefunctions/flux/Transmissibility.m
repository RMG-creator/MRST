classdef Transmissibility < StateFunction
    properties
        
    end
    
    methods
        function pp = Transmissibility(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'transMult')
                pp = pp.dependsOn({'Pressure'});
            end
        end
        
        function T = evaluateOnDomain(prop, model, state)
            T = model.operators.T;
            if isfield(model.fluid, 'transMult')
                p = prop.getEvaluatedDependencies(state, 'Pressure');
                p = model.operators.faceAvg(p);
                T = model.fluid.transMult(p).*T;
            end
        end
    end
end