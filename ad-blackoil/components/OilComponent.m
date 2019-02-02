classdef OilComponent < ImmiscibleComponent
    properties
    end
    
    methods
        function c = OilComponent(name, gasIndex)
            c@ImmiscibleComponent(name, gasIndex);
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ComponentImplementation(component, model, state);
            if model.disgas || model.vapoil
                phasenames = model.getPhaseNames();
                gix = phasenames == 'G';
                oix = phasenames == 'O';

                b = model.getProps(state, 'ShrinkageFactors');
                rhoS = model.getSurfaceDensities();
                
                if model.disgas
                    rhoOS = rhoS(oix);
                    bO = b{oix};
                    c{oix} = rhoOS.*bO;
                end
                
                if model.vapoil
                    bG = b{gix};
                    rv = model.getProp(state, 'rv');
                    c{gix} = rv.*rhoOS.*bG;
                end
            end
        end
    end
end