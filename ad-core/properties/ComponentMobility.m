classdef ComponentMobility < GridProperty
    properties

    end
    
    methods
        function v = evaluateOnDomain(prop, model, state)
            ncomp = numel(model.Components);
            v = cell(1, ncomp);
            for i = 1:ncomp
                v{i} = model.Components{i}.getComponentMobility(model, state);
            end
        end
    end
end