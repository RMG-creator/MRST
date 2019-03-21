classdef BlackOilShrinkageFactors < GridProperty
    properties
        useSaturatedFlag = false;
    end
    
    methods
        function gp = BlackOilShrinkageFactors(varargin)
            gp@GridProperty(varargin{:});
        end
        function b = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            b = cell(1, nph);
            
            f = model.fluid;
            [p, p_phase] = model.getProps(state, 'pressure', 'phasepressures');
            
            
            if model.water
                wix = phInd == 1;
                pw = p_phase{wix};
                bW = prop.evaluateFunctionOnGrid(f.bW, pw);
                b{wix} = bW;
            end
            
            if model.oil
                oix = phInd == 2;
                po = p_phase{oix};
                if model.disgas
                    rs = model.getProp(state, 'rs');
                    if prop.useSaturatedFlag
                        sG = model.getProp(state, 'sg');
                        flag = sG > 0;
                    else
                        flag = false(numelValue(p), 1);
                    end
                    bO = prop.evaluateFunctionOnGrid(f.bO, po, rs, flag);
                else
                    bO = prop.evaluateFunctionOnGrid(f.bO, po);
                end
                b{oix} = bO;
            end
            
            if model.gas
                gix = phInd == 3;
                pg = p_phase{gix};
                if model.vapoil
                    rv = model.getProp(state, 'rv');
                    if prop.useSaturatedFlag
                        sO = model.getProp(state, 'so');
                        flag = sO > 0;
                    else
                        flag = false(numelValue(p), 1);
                    end
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg, rv, flag);
                else
                    bG = prop.evaluateFunctionOnGrid(f.bG, pg);
                end
                b{gix} = bG;
            end
        end
    end
end