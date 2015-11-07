function compare_polymer(fnOPM, fnEclipse)

    mrstModule add ad-core ad-blackoil ad-fi deckformat;

    %% check/load MRST results.
    if ~(exist('resMRSTPolymer.mat', 'file') == 2)
        threePhaseBlackOilPolymerExample;
    end
    load resMRSTPolymer;

    %% check/load OPM results.
    hasOPM = (exist([fnOPM '.SMSPEC'], 'file') == 2);
    if hasOPM
        smry_opm = readSummaryLocal(fnOPM);
    end

    %% check/load MRST results.
    hasEclipse = (exist([fnEclipse '.SMSPEC'], 'file') == 2);
    if hasEclipse
        smry_ecl = readSummaryLocal(fnEclipse);
    end

    %% If we do not have results from either OPM or Eclipse, then no point to compare the results.
    if ~(hasEclipse || hasOPM)
        sprintf('There is not result from either OPM or Eclipse, the comparison will stop here!');
        return;
    end
    %%

    mrstPlotInterval = 1;

    inj = find([wellSolsPolymer{1}.sign] == 1);
    prod = find([wellSolsPolymer{1}.sign] == -1);

    injNms = {wellSolsPolymer{1}(inj).name};
    prodNms = {wellSolsPolymer{1}(prod).name};



    [qWs, qOs, qGs, bhp] = wellSolToVector(wellSolsPolymer);

    extract1 = @(wsol, fld) [ wsol.(fld) ];
    extract2 = @( c ) vertcat(c{:});
    extract  = @(fld) extract2(cellfun(@(wsol) extract1(wsol, fld), ...
                               wellSolsPolymer, 'UniformOutput', false));

    sgn =        extract('sign');
    qWPoly = sgn .* extract('qWPoly');

    Tmrst = convertTo(cumsum(schedule.step.val), day);
    if hasOPM
        Topm = smry_opm.get(':+:+:+:+', 'TIME', ':');
    end

    if hasEclipse
        Tecl = smry_ecl.get(':+:+:+:+', 'TIME', ':');
    end



    %% Polot Oil production rate
    h = figure(1);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    % set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
    mrst = qOs(:,prod);

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(prodNms{k}, 'WOPR', ':');
            plot(Topm, opm, '-g','LineWidth',2);

        end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WOPR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('m^3/day');
    title('Oil production rate');
    % axis([0 11000 0 700])

    %% Plot Water production rate
    h = figure(2);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qWs(:,prod);

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(prodNms{k}, 'WWPR', ':' );
            plot(Topm, opm, '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WWPR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('m^3/day');
    title('Water production rate');

    %% Plot Gas production rate
    h = figure(3);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qGs(:,prod);

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(prodNms{k}, 'WGPR', ':');
            plot(Topm, opm, '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WGPR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('m^3/day');
%     title('Gas production rate');

    %% Plot Water Injection rate
    h = figure(4);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qWs(:,inj);

    for k = 1:numel(inj)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(injNms{k}, 'WWIR', ':');
            plot(Topm, opm, '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(injNms{k}, 'WWIR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('m^3/day');
    title('Water injection rate');
    % axis([0 11000 0 1500]);


    %% Plot bhp for Injection well
    h = figure(5);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = bhp(:,inj);

    for k = 1:numel(inj)
        plot( Tmrst(1:mrstPlotInterval:end), convertTo(mrst(1:mrstPlotInterval:end,k), barsa), 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(injNms{k}, 'WBHP', ':');
            plot(Topm, opm, '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(injNms{k}, 'WBHP', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('bar');
    title('BHP for injection well');
    % axis([0 11000 250 450]);

    %% Plot bhp for Production well
    h = figure(6);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = bhp(:,prod);

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), convertTo(mrst(1:mrstPlotInterval:end,k), barsa), 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(prodNms{k}, 'WBHP', ':');
            plot(Topm, convertTo(opm, barsa), '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WBHP', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('bar');
%     title('BHP for production well');

    %% Plot water cut
    h = figure(7);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qWs(:,prod) ./ (qWs(:,prod) + qOs(:,prod));

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k), 'Or','LineWidth',2,'MarkerSize',13);

        if hasOPM
            opm = smry_opm.get(prodNms{k}, 'WWPR', ':') ./ ...
                 (smry_opm.get(prodNms{k}, 'WWPR', ':') + smry_opm.get(prodNms{k}, 'WOPR', ':'));
            plot(Topm, opm, '-g','LineWidth',2);
        end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WWPR', ':') ./ ...
                ( smry_ecl.get(prodNms{k}, 'WWPR', ':') + smry_ecl.get(prodNms{k}, 'WOPR', ':'));
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasOPM  && hasEclipse
        legend({'MRST', 'OPM', 'Commercial'});
    end
    if hasOPM && (~hasEclipse)
        legend({'MRST', 'OPM'});
    end
    if hasEclipse && (~hasOPM)
        legend({'MRST', 'Commercial'});
    end
    xlabel('Days');
    ylabel('Water cut');
    title('Water cut');
    % axis([0 11000 0 1]);

    %% Plot polymer for Injection well
    h = figure(8);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qWPoly(:,inj);

    for k = 1:numel(inj)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);
        % if hasOPM
        %     opm = smry_opm.get(injNms{k}, 'WCIR', ':');
        %     plot(Topm, opm, '-g','LineWidth',2);
        % end

        if hasEclipse
            ecl = smry_ecl.get(injNms{k}, 'WCIR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    % OPM can not output WCIR now.
    if hasEclipse
        legend({'MRST', 'Commercial'});
    end

    xlabel('Days');
    ylabel('kg/day');
    % axis([0 11000 0 1500]);
    title('Polymer injection rate');



    %% Plot polymer for Production well
    h = figure(9);
    set( h, 'Position', [100, 100, 900, 600]);
    clf;
    hold on;
    set(gca,'FontSize',20);
    mrst = qWPoly(:,prod);

    for k = 1:numel(prod)
        plot( Tmrst(1:mrstPlotInterval:end), mrst(1:mrstPlotInterval:end,k)*day, 'Or','LineWidth',2,'MarkerSize',13);

%         if hasOPM
%             opm = smry_opm.get(prodNms{k}, 'WCPR', ':');
%             plot(Topm, opm, '-xg');
%         end

        if hasEclipse
            ecl = smry_ecl.get(prodNms{k}, 'WCPR', ':');
            plot(Tecl, ecl, '-b','LineWidth',2,'MarkerSize',13);
        end
    end

    if hasEclipse
        legend({'MRST', 'Commercial'});
    else
        legend({'MRST'});
    end

    xlabel('Days');
    ylabel('kg/day');
    title('Polymer Production rate');

end
