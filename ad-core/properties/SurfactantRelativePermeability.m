classdef SurfactantRelativePermeability < BaseRelativePermeability

    properties
        zeroSurf
        fullSurf
    end
    
        
    methods
        function gp = SurfactantRelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            gp = gp.dependsOn({'sat', 'cs', 'capillaryNumber'}, 'state');
            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;
            gp.zeroSurf = GridProperty(model, satreg);
            gp.fullSurf = GridProperty(model, surfreg);
        end
        
        function kr = evaluateOnDomain(prop, model, state)

            fluid = model.fluid;
            [sat, cs, Nc] = model.getProps(state, 'sat', 'surfactant', 'capillaryNumber');

            isSft = (value(cs) > 0);
            m = 0*cs;
            if nnz(isSft) > 0
                logNc = log(Nc(isSft))/log(10);
                % We cap logNc (as done in Eclipse)
                logNc = min(max(-20, logNc), 20);
                m(isSft) = fluid.miscfact(logNc, 'cellInx', find(isSft));
            end
            
            [sW, sO, sG] = model.getProps(state, 'sw', 'so', 'sg');

            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;

            sWcon_noSft = fluid.krPts.w(satreg, 2);  % Residual water saturation   without surfactant
            sOres_noSft = fluid.krPts.w(satreg, 2);  % Residual oil saturation     without surfactant
            sWcon_Sft   = fluid.krPts.w(surfreg, 2); % Residual water saturation   with    surfactant
            sOres_Sft   = fluid.krPts.w(surfreg, 2); % Residual oil saturation     with    surfactant

            % Interpolated water/oil residual saturations
            sNcWcon = m.*sWcon_Sft + (1 - m).*sWcon_noSft;
            sNcOres = m.*sOres_Sft + (1 - m).*sOres_noSft;

            sNcWEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);
            sNcOEff = (sO - sNcOres)./(1 - sNcWcon - sNcOres);

            % Rescaling of the saturation - without surfactant
            sW_noSft  = (1 - sWcon_noSft - sOres_noSft).*sNcWEff + sWcon_noSft;
            sO_noSft  = (1 - sWcon_noSft - sOres_noSft).*sNcOEff + sOres_noSft;
            % Compute rel perm - without surfactant
            krW_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                              sW_noSft);
            krOW_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                              sO_noSft);
            krOG_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOG, ...
                                                              sO_noSft);
            
            d  = (sG + sW_noSft - sWcon_noSft);
            ww = (sW_noSft - sWcon_noSft)./d;
            wg = 1 - ww;
            krO_noSft = wg.*krOG_noSft + ww.*krOG_noSft;
            
            % Rescaling of the saturation - with surfactant
            sW_Sft  = (1 - sWcon_Sft - sOres_Sft).*sNcWEff + sWcon_Sft;
            sO_Sft  = (1 - sWcon_Sft - sOres_Sft).*sNcOEff + sOres_Sft;
            % Compute rel perm - with surfactant
            krW_Sft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                              sW_Sft);
            krOW_Sft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                              sO_Sft);
            krOG_Sft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOG, ...
                                                              sO_Sft);
            
            d  = (sG + sW_Sft - sWcon_Sft);
            ww = (sW_Sft - sWcon_Sft)./d;
            wg = 1 - ww;
            krO_Sft = wg.*krOG_Sft + ww.*krOG_Sft;

            % Interpolate relperm, with and without surfactant
            
            krW  = m.*krW_Sft + (1 - m).*krW_noSft;
            krO  = m.*krO_Sft + (1 - m).*krO_noSft;
            
            % We keep
            krG = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krG, sG);
            
            kr = {krW, krO, krG};
        end
    end
end