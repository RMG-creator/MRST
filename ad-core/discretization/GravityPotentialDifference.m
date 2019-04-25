classdef GravityPotentialDifference < GridProperty
    properties

    end
    
    methods
        function gp = GravityPotentialDifference(varargin)
            gp@GridProperty(varargin{:});
            gp = gp.dependsOn('Density', 'FlowPropertyFunctions');
        end
        function gRhoDz = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            gdz = model.getGravityGradient();
            if norm(model.gravity) > 0
                rho = model.getProp(state, 'Density');
                for i = 1:nph
                    rhof = model.operators.faceAvg(rho{i});
                    gRhoDz{i} = - rhof.*gdz;
                end
            else
                [gRhoDz{:}] = deal(gdz);
            end
        end
    end
end