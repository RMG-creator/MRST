classdef SurfactantRelativePermeability < BaseRelativePermeability

    properties
        zeroSurf
        fullSurf
    end
    
        
    methods
        function gp = SurfactantRelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            gp = gp.dependsOn({'sat', 'c', 'capillaryNumber'}, 'state');
            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;
            gp.zeroSurf = GridProperty(model, satreg);
            gp.fullSurf = GridProperty(model, surfreg);
        end
        
        function kr = evaluateOnDomain(prop, model, state)

            fluid = model.fluid;
            [sat, c, Nc] = model.getProps(state, 'sat', 'surfactant', 'capillaryNumber');

            isSft = (value(c) > 0);
            m = 0*c;
            if nnz(isSft) > 0
                logNc = log(Nc(isSft))/log(10);
                % We cap logNc (as done in Eclipse)
                logNc = min(max(-20, logNc), 20);
                m(isSft) = fluid.miscfact(logNc, 'cellInx', find(isSft));
            end
            
            [sW, sO, sG] = model.getProps(state, 'sw', 'so', 'sg');

            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;

            sWcon    = fluid.krPts.w(satreg, 2);  % Residual water saturation   without surfactant
            sOres    = fluid.krPts.w(satreg, 2);  % Residual oil saturation     without surfactant
            sWconSft = fluid.krPts.w(surfreg, 2); % Residual water saturation   with    surfactant
            sOresSft = fluid.krPts.w(surfreg, 2); % Residual oil saturation     with    surfactant

            % Interpolated water/oil residual saturations
            sNcWcon = m.*sWconSft + (1 - m).*sWcon;
            sNcOres = m.*sOresSft + (1 - m).*sOres;

            sNcWEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);
            sNcOEff = (sO - sNcOres)./(1 - sNcWcon - sNcOres);

            % Rescaling of the saturation - without surfactant
            sNcWnoSft  = (1 - sWcon - sOres).*sNcWEff + sWcon;
            sNcOnoSft  = (1 - sWcon - sOres).*sNcOEff + sOres;
            krNcWnoSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                              sNcWnoSft);
            krNcOnoSft = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                              sNcOnoSft);

            % Rescaling of the saturation - with surfactant
            sNcWSft  = (1 - sWconSft - sOresSft).*sNcWEff + sWconSft;
            sNcOSft  = (1 - sWconSft - sOresSft).*sNcOEff + sOresSft;
            krNcWSft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krW, ...
                                                             sNcWSft);
            krNcOSft = prop.fullSurf.evaluateFunctionOnGrid(fluid.krOW, ...
                                                             sNcOSft);

            d  = (sG + sW - sWcon);
            ww = (sW - sWcon)./d;
            wg = 1 - ww;
            
            krW  = m.*krNcWSft + (1 - m).*krNcWnoSft;

            krOW = m.*krNcOSft + (1 - m).*krNcOnoSft;
            krOG = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krOG, sO);
            krO  = wg.*krOG + ww.*krOW;
            
            krG = prop.zeroSurf.evaluateFunctionOnGrid(fluid.krG, sG);
            
            kr = {krW, krO, krG};
        end
    end
end