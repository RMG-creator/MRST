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
            sOWres_noSft = fluid.krPts.ow(satreg, 2);  % Residual oil saturation     without surfactant
            sOGres_noSft   = fluid.krPts.og(surfreg, 2); % Residual oil saturation     with    surfactant
            sGres_noSft = fluid.krPts.g(satreg, 2);  % Residual oil saturation     without surfactant
            
            sWcon_Sft   = fluid.krPts.w(surfreg, 2); % Residual water saturation   with    surfactant
            sOWres_Sft   = fluid.krPts.ow(surfreg, 2); % Residual oil saturation     with    surfactant
            sOGres_Sft   = fluid.krPts.og(surfreg, 2); % Residual oil saturation     with    surfactant
            sGres_Sft   = fluid.krPts.g(surfreg, 2); % Residual oil saturation     with    surfactant

            % Interpolated water/oil residual saturations
            sNcWcon = m.*sWcon_Sft + (1 - m).*sWcon_noSft;
            sNcOWres = m.*sOWres_Sft + (1 - m).*sOWres_noSft;
            sNcOGres = m.*sOGres_Sft + (1 - m).*sOGres_noSft;
            sNcGres = m.*sGres_Sft + (1 - m).*sGres_noSft;

            sNcWEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOWres);
            sNcOWEff = (sO - sNcOWres)./(1 - sNcWcon - sNcOWres);
            sNcGEff = (sG - sNcGres)./(1 - sNcGres - sNcOGres);
            sNcOGEff = (sO - sNcOGres)./(1 - sNcGres - sNcOGres);

            % Rescaling of the saturation - without surfactant
            sW_noSft  = (1 - sWcon_noSft - sOWres_noSft).*sNcWEff + sWcon_noSft;
            sOW_noSft  = (1 - sWcon_noSft - sOWres_noSft).*sNcOWEff + sOWres_noSft;
            sG_noSft  = (1 - sGres_noSft - sOGres_noSft).*sNcGEff + sGres_noSft;
            sOG_noSft  = (1 - sGres_noSft - sOGres_noSft).*sNcOGEff + sOGres_noSft;
            % Compute rel perm - without surfactant
            krW_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                              sW_noSft);
            krOW_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                              sOW_noSft);
            krOG_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOG, ...
                                                              sOG_noSft);
            krG_noSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krG, ...
                                                              sG_noSft);
            
            d  = (sG_noSft - sGres_noSft + sW_noSft - sWcon_noSft);
            ww = (sW_noSft - sWcon_noSft)./d;
            wg = 1 - ww;
            krO_noSft = wg.*krOG_noSft + ww.*krOW_noSft;
            
            % Rescaling of the saturation - with surfactant
            sW_Sft  = (1 - sWcon_Sft - sOWres_Sft).*sNcWEff + sWcon_Sft;
            sOW_Sft  = (1 - sWcon_Sft - sOWres_Sft).*sNcOWEff + sOWres_Sft;
            sOG_Sft  = (1 - sGres_Sft - sOGres_Sft).*sNcOGEff + sOGres_Sft;
            sG_Sft  = (1 - sGres_Sft - sOGres_Sft).*sNcGEff + sGres_Sft;
            % Compute rel perm - with surfactant
            krW_Sft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                              sW_Sft);
            krOW_Sft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                              sOW_Sft);
            krOG_Sft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krOG, ...
                                                              sOG_Sft);
            krG_Sft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krG, ...
                                                              sG_Sft);
            
            d  = (sG_Sft - sGres_Sft + sW_Sft - sWcon_Sft);
            ww = (sW_Sft - sWcon_Sft)./d;
            wg = 1 - ww;
            krO_Sft = wg.*krOG_Sft + ww.*krOW_Sft;

            % Interpolate relperm, with and without surfactant
            
            krW  = m.*krW_Sft + (1 - m).*krW_noSft;
            krO  = m.*krO_Sft + (1 - m).*krO_noSft;
            krG  = m.*krG_Sft + (1 - m).*krG_noSft;
            
            % We keep
%             krG = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krG, sG);
            
            kr = {krW, krO, krG};
        end
    end
end