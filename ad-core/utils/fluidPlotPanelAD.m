function fluidPlotPanelAD(model, varargin)
% Create simple plot of fluid model for a given simulation model.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    opt = struct('pressureRange', []);
    opt = merge_options(opt, varargin{:});

    if isempty(opt.pressureRange)
        p0 = max(model.minimumPressure, 0);
        p1 = min(model.maximumPressure, 1000*barsa);
        opt.pressureRange = subdiv(p0, p1);
    end
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    fh = figure('Position', df.*[1 1 1.75 1]);
    % We want to ensure that the lines are nice and pretty even if the user
    % has messed with the defaults.
    set(fh, 'Renderer', 'painters');
    
    % Options can alter the amount of space the plot itself takes up.
    lm = 0.1;
    pw = .7;
    % Somewhat magic numbers because the gui in matlab has some magic
    % constants itself.
    plotaxis  = subplot('Position', [.75*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, .25*lm, 1-1.25*lm-pw,  1-1.25*lm], ...
                        'Title',  'Property');
    
    names = {'Viscosity',...
             'Relperm', ...
             'Rock compressibility', ...
             'Densities', ...
             'Capillary pressure', ...
             'Max Rs/Rv'
             };
    functions = {@(name) plotViscosity(model, name), ...
                 @(name) plotRelperm(model, name),...
                 @(name) plotPVMult(model, name),...
                 @(name) plotDensity(model, name), ...
                 @(name) plotCapillaryPressure(model, name), ...
                 @(name) plotMaxR(model, name)
                };
    
    if model.oil
        names = [names, 'Oil viscosity'];
        names = [names, 'Oil b-factor'];
    end
    
    
%     if isfield(model.fluid, 'pcOW') || isfield(model.fluid, 'pcOG')
%         names{end + 1} = 'Capillary pressure';
%     end
%     if checkBO(model)
%         names{end + 1} = 'Max dissolution';
%     end
    propsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox',...
              'String', names, 'Callback', @drawPlot, ...
              'Position',[0 0 1 1]);
    
          

    function drawPlot(src, event)
        axis(plotaxis);
        ix = get(propsel, 'Value');
        fn = functions{ix};
        name = names{ix};
        
        fn(name);
    end

    drawPlot([], []);


    function plotStuff(model, fields)
        f = model.fluid;
        p = opt.pressureRange;
        s = subdiv(0, 1);

        rs = 0;
        rv = 0;
        if model.disgas
            rs = f.rsSat(p);
        end
        if model.vapoil
            rv = f.rvSat(p);
        end
        n = sum(model.getActivePhases);

        legflag = false(size(fields));
        data = nan(numel(p), n);

        ctr = 0;
        yl = '';
        for i = 1:numel(fields)
            fn = fields{i};
%             if ~isfield(f, fn)
%                 continue
%             end
            legflag(i) = true;
            ctr = ctr + 1;
            switch(lower(fn))

                case {'krw', 'krg', 'krog', 'kro', 'krow'}
                    data(:, ctr) = f.(fn)(s);
                    x = s;
                    xl = 'Saturation';
                case {'pcow', 'pcog'}
                    data(:, ctr) = f.(fn)(s);
                    x = s;
                    xl = 'Saturation';
                case {'bo', 'bg', 'bw'}
                    data(:, ctr) = evalSat(model, f, fn, p, rs, rv);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
                case {'rhoo', 'rhog', 'rhow'}
                    bsub = ['b', fn(end)];
                    rho = f.(['rho', fn(end), 'S']);
                    b = evalSat(model, f, bsub, p, rs, rv);
                   
                    data(:, ctr) = b*rho;
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
                    yl = 'Density [kg/m^3]';
                case {'muw', 'muo', 'mug'}
                    data(:, ctr) = evalSat(model, f, fn, p, rs, rv);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
                case {'rssat', 'rvsat'}
                    data(:, ctr) = f.(fn)(p);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
                case {'pvmultr'}
                    data(:, ctr) = f.(fn)(p);
                    x = p/barsa;
                    xl = 'Pressure (barsa)';
            end
        end
        plot(x, data, 'linewidth', 2)
        grid on
        legend(fields(legflag))
        xlabel(xl)
        ylabel(yl);
    end

    function plotViscosity(model, name)
        plotStuff(model, {'muW', 'muO', 'muG'});
        if checkBO(model)
            title([name, ' (saturated curves)']);
        else
            title([name, 'Viscosity']);
        end
    end

    function plotRelperm(model, name)
        krnames = {};
        if model.water
            krnames = [krnames, 'krW'];
        end
        if model.oil
            if model.water && model.gas
                krnames = [krnames, 'krOW', 'krOG'];
            else
                krnames = [krnames, 'krO'];
            end
        end
        if model.gas
            krnames = [krnames, 'krG'];
        end
        plotStuff(model, krnames);
        title(name);
    end

%     function plotCompressibility(model, name)
%         plotStuff(model, {'bW', 'bO', 'bG'});
%         if checkBO(model)
%             title([name, ' (saturated curves)']);
%         else
%             title(name);
%         end
%     end
    function plotDensity(model, name)
        plotStuff(model, {'rhoW', 'rhoO', 'rhoG'});
        if checkBO(model)
            title([name, ' (saturated curves)']);
        else
            title(name);
        end
    end

    function plotCapillaryPressure(model, name)
        plotStuff(model, {'pcOW', 'pcOG'});
        title('Capillary pressure');
    end

    function plotMaxR(model, name)
        rnames = {};
        if model.disgas
            rnames = [rnames, 'rsSat'];
        end
        if model.vapoil
            rnames = [rnames, 'rvSat'];
        end
        plotStuff(model, rnames);
        title(name);
    end

    function plotPVMult(model, name)
        plotStuff(model, {'pvMultR'});
        title(name);
    end

    function x = subdiv(start, stop)
        dx = (stop - start)/250;
        x = (start:dx:stop)';
    end

    function y = evalSat(model, f, fn, x, rs, rv)
        if checkBO(model)
            if any(strcmpi(fn, {'muo', 'bo'})) && model.disgas
                y = f.(fn)(x, rs, true(size(x)));
            elseif any(strcmpi(fn, {'mug', 'bg'})) && model.vapoil
                y = f.(fn)(x, rv, true(size(x)));
            else
                y = f.(fn)(x);
            end
        else
            y = f.(fn)(x);
        end
    end

    function ind = checkBO(model)
        ind = isa(model, 'ThreePhaseBlackOilModel') &&...
               (model.disgas || model.vapoil);
    end
end