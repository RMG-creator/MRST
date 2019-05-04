classdef FaceConcentration < GridProperty & UpwindProperty
    properties
        
    end
    
    methods
        function gp = FaceConcentration(backend, upwinding)
            gp@GridProperty(backend);
            gp@UpwindProperty(upwinding)
            gp = gp.dependsOn('PhaseUpwindFlag');
        end
        
        function fc = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, 'PhaseUpwindFlag');
            c = model.getProps(state, 'polymer');
            fc = prop.faceUpstream(state, flag{1}, c);
            
        end
    end
end