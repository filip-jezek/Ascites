%% Load paracentesis
dataTable = readtable('../data/pressure-volume-time_Adj.csv')

plot(dataTable.Var1, dataTable.Liters_Drained__, dataTable.Var1,dataTable.Record_ID_Adj);
%% Process data into patient struct
patient = struct();
patients = {};
rec = struct();
for i_row = 1:height(dataTable)
    line = dataTable(i_row, :);
    patId = line.Subject_ID_Adj;
    recId = line.Record_ID_Adj;
    rec.datedate = line.Date_of_paracentesis_;
    
    if length(patients) < patId || isempty(patients{patId})
        patients{patId} = struct();
    end
    if ~isfield(patients{patId}, 'rec') || length(patients{patId}.rec) < recId
        patients{patId}.rec{recId} = struct();
    end
    if ~isfield(patients{patId}.rec{recId}, 'Drained')
        patients{patId}.rec{recId}.Drained = [];
        patients{patId}.rec{recId}.Pressure = [];
    end
    patients{patId}.rec{recId}.Drained = [patients{patId}.rec{recId}.Drained line.Liters_Drained__];
    patients{patId}.rec{recId}.Pressure = [patients{patId}.rec{recId}.Pressure line.Pressure_Measurement__];
    patients{patId}.rec{recId}.Date = line.Date_of_paracentesis_;
    patients{patId}.rec{recId}.last2dop = line.last2dop;
    if isnan(line.last2dop)
        patients{patId}.FirstParaEver = true;
    end
end
% squish empty procedure    
isSame = 0;isLarger = 0;isLower = 0;isNan = 0;
for patId = 1:length(patients)
    if isempty(patients{patId})
        continue;
    end
    for recId = 1:length(patients{patId}.rec)
        if isempty(patients{patId}.rec{recId})
            continue;
        end
        if isfield(patients{patId}, 'prc')
            prcPos = length(patients{patId}.prc) + 1;
        else
            prcPos = 1;
        end
        patients{patId}.prc{prcPos} = patients{patId}.rec{recId};
        patients{patId}.dayOfPrc(prcPos) = patients{patId}.rec{recId}.Date;
        patients{patId}.last2dop(prcPos) = patients{patId}.rec{recId}.last2dop;
        if isnan(patients{patId}.last2dop(prcPos))
            isNan = isNan + 1;
        elseif prcPos > 1 
            % check the data
            if patients{patId}.rec{recId}.last2dop > days(patients{patId}.dayOfPrc(prcPos) - patients{patId}.dayOfPrc(prcPos-1))
                isLarger = isLarger + 1;
                % correct the data
                patients{patId}.last2dop(prcPos) = days(patients{patId}.dayOfPrc(prcPos) - patients{patId}.dayOfPrc(prcPos-1));
            elseif patients{patId}.rec{recId}.last2dop < days(patients{patId}.dayOfPrc(prcPos) - patients{patId}.dayOfPrc(prcPos-1))
                isLower = isLower + 1;
            else
                isSame = isSame + 1;
            end
        end
    end    
    patients{patId} = rmfield(patients{patId}, "rec");
end
sprintf('It is same for %d procedures, \nshorter for %d (missing proc) and \nlonger (lying) for %d procedures', isSame, isLower, isLarger)
%% plot that
figure(1);clf;
tiledlayout('flow');

for patId = 1:length(patients)
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
%%
    nexttile;
    
    for prcId = 1:length(patients{patId}.prc)
        prc = patients{patId}.prc{prcId};
        plot(prc.Drained, prc.Pressure, 'o-');hold on;
        
    end
    title(sprintf('Patient %d', patId));
    xlim([0 inf]);ylim([0 inf]);
end
%%
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.01, 'MaxIter', 20);
% 
% %% single patient MODEL fit and estimation
% figure(2);clf;
% tiledlayout('flow');
% nexttile;
% 
% % linear model
% pm = @(a, b, c, x) max(-1, a*(x -b) + c);
% % no shift
% pm_0 = @(a, c, x) pm(a, 0, c, x);
% costs = [];
% % patId = 2;
% % patId = 2;
% patId = 11;
% % prcId = 1;
% % patId = 10;
% x_ = 0:10;
% hold on;
% cl = lines(length(patients{patId}.prc));
% for prcId = 1:size(patients{patId}.prc, 2)
%     if prcId == 1
%         lw = 3; % line width
%     else
%         lw = 1;
%     end
%     % if patients{patId}.prc{prcId}.Pressure == 0
%     %     patients{patId}.prc(prcId).Pressure = NaN;
%     % end
%     % [ae gd] = fit(patients{patId}.prc{prcId}.Drained', patients{patId}.prc{prcId}.Pressure', pm, ...
%     %     'StartPoint', [-1, 0, patients{patId}.prc{prcId}.Pressure(1)]);
%     [ae gd] = fit(patients{patId}.prc{prcId}.Drained', patients{patId}.prc{prcId}.Pressure', pm_0, ...
%         'StartPoint', [-1, patients{patId}.prc{prcId}.Pressure(1)]);
% 
%     costs(prcId) = gd.rmse/length(patients{patId}.prc{prcId}.Drained);
%     fsoptions = optimset('Display','off');
%     patients{patId}.Vdrained(prcId) = patients{patId}.prc{prcId}.Drained(end) - patients{patId}.prc{prcId}.Drained(1);
%     V0 = fsolve(ae, 1, fsoptions);
%     x_ = 0:0.1:ceil(V0);
%     patients{patId}.V0(prcId) = V0;
%     patients{patId}.V_openExtrap(prcId) = V0 - patients{patId}.prc{prcId}.Drained(1);
%     patients{patId}.V_closeExtrap(prcId) = V0 - patients{patId}.prc{prcId}.Drained(end);
%     patients{patId}.V_closeExtrapByP(prcId) = V0 - fsolve(@(x) ae(x) - patients{patId}.prc{prcId}.Pressure(end), 1, fsoptions);
%     if prcId == 1
%         patients{patId}.V_generated(prcId) = NaN;
%     else
%         patients{patId}.V_generated(prcId) = patients{patId}.V_openExtrap(prcId) - patients{patId}.V_closeExtrap(prcId - 1);
%     end
% 
%     plot(V0 - patients{patId}.prc{prcId}.Drained, patients{patId}.prc{prcId}.Pressure, 'o-', 'Color', cl(prcId, :), LineWidth=lw);
%     plot(V0 - x_, ae(x_), '--', LineWidth=lw, Color=cl(prcId, :));
% end
% cost = nanmean(costs);
% 
% title(sprintf('Single procedure best fit: %0.2fx + %0.2f', ae.a, ae.c))
% xlabel('$V_E = V_0 - V_D$', 'Interpreter','latex');ylabel('Pressure ($cmH_{2}O$)', 'Interpreter','latex');
% %% Fit all procedures in single patient
% nexttile;hold on;
% % fitPat = @(params) evalFit(params, patients{patId}.prc, false);
% % fr = fminsearch(fitPat, [-1 10], options);
% %% fit only zero crossings
% % clf;
% fitPat = @(params) evalFit(params, patients{patId}.prc, false);
% fr = fminsearch(fitPat, [10], options);
% 
% % cla;hold on;
% [costAvg bs costs] = evalFit([fr], patients{patId}.prc, true);
% costAvg
% %
% if length(fr) == 1
%     title(sprintf('All procedures best fit: V0 at %0.2f L', fr(1)));
% else
%     title(sprintf('All procedures best fit: %0.2f(x - V_0) + %0.2f', fr(1), fr(2)));
% end
% xlabel('$V_D - V_{S,j}$', 'Interpreter','latex');ylabel('Pressure ($cmH_2O$)', 'Interpreter','latex');
% %% Visualize date shift
% max drained volume, which is empty ambdomen
% pm_zc
% V_dmax = fsolve(@(x)pm(fr(1), 0, fr(2), x), 1e-6);
% 
% for prcId = 1:length(bs)
%     prc = patients{patId}.prc{prcId};
%     patients{patId}.V_openCombined(prcId) = -(prc.Drained(1) - V_dmax - bs(prcId));
%     plot(prc.Drained - V_dmax - bs(prcId), prc.Pressure, 'o-', Color=cl(prcId, :));hold on;
% end
% x_ = 0:0.1:10;
% plot(x_ - V_dmax, pm(fr(1), 0, fr(2), x_))
% plot([V_dmax V_dmax], [0  10], 'k--');
% text(V_dmax, 10, sprintf('$V_{D,max}$\n %0.1f L',V_dmax), 'Interpreter','latex')
% patients{patId}.V_dmax = V_dmax;
% days from last paracentesis
% dflp = patients{patId}.last2dop;
% nexttile;hold on;
% plotXXXonDFLP(dflp, patients{patId}.Vdrained, "Drained volume", []);
% nexttile;hold on;
% plotXXXonDFLP(dflp, patients{patId}.V_openExtrap, "Opening volume - Extrapolation", [])
% nexttile;hold on;
% plotXXXonDFLP(dflp, patients{patId}.V_closeExtrap, "Closing volume - extrapolation", [])
% nexttile;hold on;
% plotXXXonDFLP(dflp, patients{patId}.V_generated, "Generated volume", [])
% 



%% Valid patients - shrink total size
allPatIds = 1:min(120,length(patients));
% allPatIds = [2 3];
patIds = [];
for patId = allPatIds
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
    patIds = [patIds, patId];
end

%% valid procedures and outliers
for patId = patIds
    validPrcs = false(size(patients{patId}.prc));
    for prcId = 1:length(patients{patId}.prc)
        % datapoint is valid, if the pressure is monotonically increasing
        % from zero. Any decrease is discarded
        % flipping it to diff upwards and then back
        patients{patId}.prc{prcId}.validDP = fliplr([true diff(fliplr(patients{patId}.prc{prcId}.Pressure)) > 0]);
        len = sum(patients{patId}.prc{prcId}.validDP);

        if len >= 2
            % we need at least two points to identify the line
            validPrcs(prcId) = true;
        else
            validPrcs(prcId) = false;
        end
    end
    patients{patId}.prc = patients{patId}.prc(validPrcs);
    patients{patId}.last2dop = patients{patId}.last2dop(validPrcs);
end

%% Identify all patients - assuming same compliance
reoptim = true;
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.01, 'MaxIter', 10);
fsoptions = optimset('Display','off');

figure(3);clf;
params = zeros(length(patId), 2);
initParams = [-1 10];

for patId = patIds(1:end) % patIds
    nexttile;hold on;
    prcs = patients{patId}.prc;
    fitPat = @(params) evalFit(params, prcs, false);
    
    if reoptim
        fr = fminsearch(fitPat, initParams, options);
        [costAvg, shift, costs] = evalFit(fr, prcs, true);
        patients{patId}.params = fr;
        patients{patId}.shifts = shift;
        patients{patId}.costs = costs;
    else
        [costAvg, shift, costs] = evalFit(patients{patId}.params, prcs, true);
        fr = patients{patId}.params;
    end

    
    % max drained volume, which is zero-pressure volume
    % we get it as a drained volume, at which the model crosses 0 pressure
    % (from above, thus smallest positive number will do)
    
    c = lines(length(patients{patId}.prc));
    for prcId = 1:min(inf, length(prcs))
        % if length(patients{patId}.shifts) < prcId
        %     % last procedure, unable to fit sinlge data point
        %     continue;
        % end
        V0 = patients{patId}.shifts(prcId);
        
        patients{patId}.v0Comb(prcId) = V0;
        prc = patients{patId}.prc{prcId};
        p_open = prc.Pressure(1);
        patients{patId}.p_open(prcId) = p_open;                

        plot(V0 - prc.Drained, prc.Pressure, 's--')

        % vest_at_open = fsolve(@(x)pm(fr(1), 0, fr(2), x) - p_open, 1e-6, fsoptions);
        
        patients{patId}.Vdrained(prcId) = prc.Drained(end) - prc.Drained(1);

        patients{patId}.V_openCombined(prcId) = V0 - prc.Drained(1);        
        patients{patId}.V_closeCombined(prcId) = V0 - prc.Drained(end);
        
        if prcId == 1
            patients{patId}.V_generatedComb(prcId) = NaN;
        else
            patients{patId}.V_generatedComb(prcId) = patients{patId}.V_openCombined(prcId) - patients{patId}.V_closeCombined(prcId - 1);
            plot(patients{patId}.V_closeCombined(prcId-1), patients{patId}.prc{prcId-1}.Pressure(end), '>', MarkerSize=12, Color=c(prcId, :), LineWidth=3)
            plot(patients{patId}.V_openCombined(prcId), prc.Pressure(1), '<', MarkerSize=12, Color=c(prcId, :), LineWidth=3)           
        end
    end

    % params(patId, :) = fr;
    title(sprintf('Patient %d costs %0.2f', patId, costAvg));
    xlim([-inf inf]);ylim([0 inf]);    
end

%% METHOD 2: Extrapolating for zero pressure
pm = @(a, b, x) max(0, a*(x - b));
figure(4);clf; tiledlayout('flow');
for patId = patIds

    nexttile;hold on;
    fsoptions = optimset('Display','off');
    prc = patients{patId}.prc;
    fitPat = @(params) evalFit(params, prc, false);
    
    c = lines(length(patients{patId}.prc));
    for prcId = 1:min(inf, length(patients{patId}.prc))
        prc = patients{patId}.prc{prcId};
        if length(patients{patId}.prc{prcId}.Drained) < 2
            continue;
        end
        
        if sum(prc.validDP) < 2
            % disregeard if less than 2 valid data points
            continue;
        end
        % identify by extrapolating each one         
        [ae gd] = fit(prc.Drained(prc.validDP)', prc.Pressure(prc.validDP)', pm, ...
            'StartPoint', [-1, prc.Pressure(1)]);        
        costs(prcId) = gd.rmse/(length(prc.Drained(prc.validDP)));
        patients{patId}.costsExtrap(prcId) = costs(prcId);
        V0 = fsolve(@(x)ae(x) - 1e-6, 1, fsoptions);
        x_ = 0:0.1:ceil(V0);
        patients{patId}.V0(prcId) = V0;
        patients{patId}.V_openExtrap(prcId) = V0 - prc.Drained(1);
        patients{patId}.ExtrapSlope(prcId) = ae.a;

        patients{patId}.V_closeExtrap(prcId) = V0 - prc.Drained(end);
        patients{patId}.V_closeExtrapByP(prcId) = V0 - fsolve(@(x) ae(x) - prc.Pressure(end), 1, fsoptions);
        if prcId == 1
            patients{patId}.V_generatedExtrap(prcId) = NaN;
        else
            patients{patId}.V_generatedExtrap(prcId) = patients{patId}.V_openExtrap(prcId) - patients{patId}.V_closeExtrap(prcId - 1);
            plot(patients{patId}.V_closeExtrap(prcId-1), patients{patId}.prc{prcId-1}.Pressure(end), '>', MarkerSize=12, Color=c(prcId, :), LineWidth=3)
            plot(patients{patId}.V_openExtrap(prcId), prc.Pressure(1), '<', MarkerSize=12, Color=c(prcId, :), LineWidth=3)
        end

        plot(patients{patId}.V0(prcId) - prc.Drained(prc.validDP), prc.Pressure(prc.validDP), 'o-', Color=c(prcId, :), LineWidth=2, markerSize=8)
        plot(patients{patId}.V0(prcId) - prc.Drained(~prc.validDP), prc.Pressure(~prc.validDP), 'x', Color=c(prcId, :), LineWidth=4, markerSize=12)
        x_ = 0:0.1:20;
        plot(patients{patId}.V0(prcId) - x_, ae(x_), '--', Color=c(prcId, :));

    end

    % params(patId, :) = fr;
    title(sprintf('Patient %d costs %0.2f', patId, nanmean(costs)));
    xlim([-inf inf]);ylim([0 inf]);
end

%% PLOTS - correlation between combined and etrapolated drainage
figure(30);clf;tiledlayout('flow');
c = zeros(max(patIds), 3); c(patIds, :) = lines(length(patIds));
nexttile;
hold on;
plot([0 15], [0 15], 'k--')
for patId = patIds
    plot(patients{patId}.V_generatedExtrap, patients{patId}.V_generatedComb, 'o', LineWidth=2, Color=c(patId, :));
    pl(patId) = plot(patients{patId}.V_generatedExtrap, patients{patId}.V_generatedComb, 'o-', LineWidth=0.5, Color=c(patId, :));
    leg(patId) = string(sprintf('#%d', patId));
    valid(patId) = true;
end
xlim([0 inf]);ylim([0 inf]);
xlabel('Extrap');ylabel('Combined');legend(pl(valid), leg(valid))
%% COMBINED
figure(31);clf;tiledlayout('flow', TileSpacing='compact');

% COMBINED
nexttile;hold on;
for patId = patIds
    pl(patId) = plot(patients{patId}.V_generatedComb, 'o-', LineWidth=2);
end
ylim([0 inf]);
xlabel('para Sequence');ylabel('V_{gen} Combined');%legend(pl(valid), leg(valid));

nexttile;hold on;
for patId = patIds
    if length(patients{patId}.V_generatedComb) < 2
        continue;
    end
    ps = plot(patients{patId}.last2dop(2), patients{patId}.V_generatedComb(2), 's', LineWidth=4, color=c(patId, :));    
    pl(patId) = plot(patients{patId}.last2dop, patients{patId}.V_generatedComb, 'o-', LineWidth=2, color=c(patId, :));
end
ylim([0 inf]);
xlabel('Days since last para');ylabel('V_{gen} Combined');%legend(pl(valid), leg(valid));

nexttile;hold on;
for patId = patIds
    pl(patId) = plot(patients{patId}.V_generatedComb./patients{patId}.last2dop, 'o-', LineWidth=2);
end
ylim([0 inf]);
xlabel('para Sequence');ylabel('R_{gen} Combined (L/day)');%legend(pl(valid), leg(valid));

nexttile;hold on;
for patId = patIds
    if length(patients{patId}.V_generatedComb) < 2
        continue;
    end    
    ps(patId) = plot(patients{patId}.last2dop(2), patients{patId}.V_generatedComb(2)./patients{patId}.last2dop(2), 's', LineWidth=4, color=c(patId, :));
    pl(patId) = plot(patients{patId}.last2dop, patients{patId}.V_generatedComb./patients{patId}.last2dop, 'o-', LineWidth=2, color=c(patId, :));
end
ylim([0 inf]);
xlabel('Days since last para');ylabel('R_{gen} Combined (L/day)');legend(pl(valid), leg(valid));
exportgraphics(gcf, 'GeneratedVolumes_combined.png', 'Resolution',150);

%% EXTRPOLATED
figure(32);clf;tiledlayout('flow', 'TileSpacing','compact', 'Padding','loose');
% fill in colors only for valid patients
c = zeros(max(patIds), 3); c(patIds, :) = lines(length(patIds));
nexttile;hold on;
for patId = patIds
    pl(patId) = plot(patients{patId}.V_generatedExtrap, 'o-', LineWidth=2);
end
ylim([0 inf]);
xlabel('para Sequence');ylabel('V_{gen} Extrap');%legend(pl(valid), leg(valid));

nexttile;hold on;
for patId = patIds
    if length(patients{patId}.V_generatedExtrap) < 2
        continue;
    end
    ps = plot(patients{patId}.last2dop(2), patients{patId}.V_generatedExtrap(2), 's', LineWidth=4, color=c(patId, :));
    pl(patId) = plot(patients{patId}.last2dop, patients{patId}.V_generatedExtrap, 'o-', LineWidth=2, color=c(patId, :));
end
ylim([0 inf]);
xlabel('Days since last para');ylabel('V_{gen} Extrap');%legend([pl(valid) ps], [leg(valid) "Start"]);

nexttile;hold on;
for patId = patIds
    pl(patId) = plot(patients{patId}.V_generatedExtrap./patients{patId}.last2dop, 'o-', LineWidth=2);
end
ylim([0 inf]);
xlabel('para Sequence');ylabel('R_{gen} Extrap (L/day)');%legend(pl(valid), leg(valid));

nexttile;hold on;
for patId = patIds
    if length(patients{patId}.V_generatedExtrap) < 2
        continue;
    end
    ps(patId) = plot(patients{patId}.last2dop(2), patients{patId}.V_generatedExtrap(2)./patients{patId}.last2dop(2), 's', LineWidth=4, color=c(patId, :));
    pl(patId) = plot(patients{patId}.last2dop, patients{patId}.V_generatedExtrap./patients{patId}.last2dop, 'o-', LineWidth=2, color=c(patId, :));
end
ylim([0 inf]);
xlabel('Days since last para');ylabel('R_{gen} Extrap (L/day)');legend(pl(valid), leg(valid));
exportgraphics(gcf, 'GeneratedVolumes_extrapolated.png', 'Resolution',150);
%% PLOTS - drained volume
maxPat = inf;
figure(5);clf;tl = tiledlayout('flow');
title(tl,'Drained volume by days of last paracentesis')
for patId = 1:min(maxPat, length(patients))
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
    nexttile; hold on;
    dflp = patients{patId}.last2dop;
    plotXXXonDFLP(dflp, patients{patId}.Vdrained, "V_D", []);
    title(sprintf('Pat #%d', patId));l = legend();set(l,'visible', 'off');
end
%%
figure(6);clf;tl = tiledlayout('flow');
title(tl,'Slope of (Extrapolated) by days of last paracentesis')
for patId = 1:min(maxPat, length(patients))
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
    nexttile; hold on;
    dflp = patients{patId}.last2dop;
    plotXXXonDFLP(dflp, -patients{patId}.ExtrapSlope, "-Slope_{PV}", []);
    title(sprintf('Pat #%d', patId));l = legend();set(l,'visible', 'off');
end
%%
figure(6);clf;tl = tiledlayout('flow');
title(tl,'Generated volume')
for patId = 1:min(maxPat, length(patients))
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
    nexttile; hold on;
    dflp = patients{patId}.last2dop;
    plotXXXonDFLP(dflp, patients{patId}.V_generated, "V_G", []);
    title(sprintf('Pat #%d', patId));l = legend();set(l,'visible', 'off');
end
%% plot correlations
figure(6);clf;
tl = tiledlayout('flow');
title(tl, 'Some correlations')

ax1 = nexttile;hold on;
xlabel('Opening volume');ylabel('Slope (comb)');

ax2 = nexttile;hold on;
xlabel('Drained volume');ylabel('Slope (comb)');

ax3 = nexttile;hold on;
xlabel('Drained volume');ylabel('Slope (Extrap)');

ax4 = nexttile;hold on;
xlabel('Time since last para');ylabel('Slope (Extrap)');

ax5 = nexttile;hold on;
xlabel('Time since last para');ylabel('Openinig pressure');

ax6 = nexttile;hold on;
xlabel('Time since last para');ylabel('Openinig volume (Comb)');

ax7 = nexttile;hold on;
xlabel('Time since last para');ylabel('Openinig volume (Extrap)');

ax8 = nexttile;hold on;
xlabel('Time since last para');ylabel('Drained volume');

ax9 = nexttile;hold on;
xlabel('V_{open} (Comb)');ylabel('V_{open} (Extrap)');xlim([0, 20]);ylim([0 20])


for patId = 1:length(patients)
    if ~isfield(patients{patId}, 'params')
        continue;
    end
    plot(ax1, patients{patId}.V_openCombined , patients{patId}.params(1), 'o')
    plot(ax2, patients{patId}.Vdrained, patients{patId}.params(1), 'o')
    plot(ax3, patients{patId}.Vdrained, patients{patId}.ExtrapSlope, 'o')
    plot(ax4, patients{patId}.last2dop, patients{patId}.ExtrapSlope, 'o')
    plot(ax5, patients{patId}.last2dop, patients{patId}.p_open, 'o'); 
    plot(ax6, patients{patId}.last2dop, patients{patId}.V_openCombined, 'o'); 
    plot(ax7, patients{patId}.last2dop, patients{patId}.V_openExtrap, 'o'); 
    plot(ax8, patients{patId}.last2dop, patients{patId}.Vdrained, 'o'); 
    plot(ax9, patients{patId}.V_openCombined, patients{patId}.V_openExtrap, 'o'); 
end
% ylim([0 inf]);xlim([0 inf])
%% save to json
fid  = fopen('../data/para_output.json', 'w');
fprintf(fid, jsonencode(patients, 'PrettyPrint', true))
fclose(fid);
%% save to csv
names = ["Patient", "Procedure", ...
    "DaysSinceLastPara", "TotalDrained", "P_open", "P_close", ...
    "V_open_comb", "V_close_comb", "V_gen_comb", "slope_comb", "rate_comb","cost_comb"...
    "V_openExtrap", "V_closeExtrap", "V_generatedExtrap", "ExtrapSlope", "rate_extrap","cost_extrap"];
t = table;
for patId = patIds
    pat = patients{patId};
    for prcId = 1:length(patients{patId}.prc)
        prc = patients{patId}.prc{prcId};
        
        % first 17 patients in feasibility study
        feasMax = 17;
        if patId > feasMax
            % follow/up study
            patIdAdj = 100+patId - feasMax;
        else
            % feasibility study
            patIdAdj = patId;
        end
        tr = {patIdAdj, prcId, ...
            pat.last2dop(prcId), pat.Vdrained(prcId), pat.p_open(prcId), prc.Pressure(end), ...
            pat.V_openCombined(prcId), pat.V_closeCombined(prcId), pat.V_generatedComb(prcId), -pat.params(1), pat.V_generatedComb(prcId)./pat.last2dop(prcId), pat.costs(prcId), ...
            pat.V_openExtrap(prcId), pat.V_closeExtrap(prcId), pat.V_generatedExtrap(prcId), pat.ExtrapSlope(prcId), pat.V_generatedComb(prcId)./pat.last2dop(prcId), pat.costsExtrap(prcId)};
        t(end + 1, :) = tr;
    end
end
t.Properties.VariableNames = names;
% T = table('Size', [0 2], 'VariableTypes', {'double', 'double'}, 'VariableNames', names);
writetable(t, '../data/ExportParacentesis.csv');

%% Funciton evalFit
function [costAvg fitparam costs]= evalFit(params, prcSet, showPlots)
    v0 = 0;
    % if true first procedure is not shifted
    useReference = false;
    % if true we disregard all zero pressure values for fitting
    useZeroCutOff = false;
    
    if length(params) == 2
    % identify only shift of each procedure here
    % slope and vertical shift is common
    a = params(1);
    c = params(2);
    c = 0;
    pm_base = @(b, x) max(0, a*(x - b) + c);
    pm_0 = @(b, x) pm_base(0, x) + 0*b;
    init = [10];
    elseif length(params) == 1
        % only zero crossing is common
        zc = params(1);
        % a(x - b)+c = 0;
        % x = b -c/a;
        % zc = 9;
        % a(zc-b) + c = 0;
        % c = -a(zc-b)
        pm_base = @(a, b, x) max(0, a*(x - b) -a*zc);        
        pm_0 = @(a, b, x) pm_base(a, 0,x) + 0*b;
        init = [-1,0];
    % elseif isempty(params)
    %     % identify all, nothing in common
    %     pm_base = @(a, b, c, x) max(0, a*(x - b*0) + c);
    %     pm_0 = pm_base;
    %     init = [-1, 0, 1];
    end
    %%
    costs = zeros(length(prcSet), 1);
    fitparam = [];
    
    % reference fit - just for the costs here
    % [ae gd] = fit(prcSet{1}.Drained', prcSet{1}.Pressure', pm, 'StartPoint', [0]);
    % cost = gd.rmse;
    
    if showPlots
    %     plot(prcSet{1}.Drained, prcSet{1}.Pressure);
    % 
    %     hold on;
        x_ = 0:0.01:10;
    %     plot(x_, pm_base(0, x_), '--')
    end
    cl = lines(length(prcSet));
    %%
    for i = 1:length(prcSet)
        if useReference && i == 1
            % all shifts relative to first one
            pm = pm_0;
        else
            pm = pm_base;
        end
        prc = prcSet{i};
        
        if useZeroCutOff
            % experimental: remove 0 values
            valid = (prc.Pressure > 0);
            prc.Drained = prc.Drained(valid);
            prc.Pressure = prc.Pressure(valid);
        end

        
        if sum(prc.validDP) < 2
            % disregeard if less than 2 valid data points
            continue;
        end
        % shift other drains
        [ae gd] = fit(prc.Drained(prc.validDP)', prc.Pressure(prc.validDP)', pm, 'StartPoint', init);        
        costs(i) = gd.rmse/(length(prc.Drained(prc.validDP)));
        fitparam(i, :) = coeffvalues(ae);        
        if showPlots
            % v0 = fsolve(@(x)ae(x) - 1e-3, 1); == b :)
            v0 = ae.b;
            pleg(i) = plot(v0 - prc.Drained(prc.validDP), prc.Pressure(prc.validDP), 'o-', Color=cl(i, :), LineWidth=2, MarkerSize=8);
            hold on;
            plot(v0 - prc.Drained(~prc.validDP), prc.Pressure(~prc.validDP), 'x', Color=cl(i, :), LineWidth=4, MarkerSize=12)
            x_ = 0:0.01:10;
            plot(v0 - x_, ae(x_), '--', LineWidth=2, Color=cl(i, :));
            % leg(i) = string(sprintf('V_{0} %0.2f', fitparam(i)));
        end
    end
    %%
    costAvg = mean(costs);
    if showPlots
        % legend(pleg, leg, 'AutoUpdate','off')
    end
end
%% function plotXXX
function succ = plotXXXonDFLP(dflp, quantity, quantityTitle, ms)
if length(dflp) ~= length(quantity)
    disp('Arrays are not the same, quitting here')
    return;
end
cl = lines(length(dflp));
    if isempty(ms)
        ms = 2*ones(length(dflp));
    end
    if isnan(dflp(1))
        % first paracentesis ever
        lp(1) = plot([0 max(dflp(2:end))*1.5], [quantity(1) quantity(1)], '--', LineWidth=2, Color=cl(1, :));hold on;
        lpl(1) = string("#1 (last unknown)");
        first = 2;
    else
        first = 1;
    end
    % first paracentesis
    % cla;
    % clear lp lpl;
    plot(dflp(first:end), quantity(first:end), 'k--');
    for prcId = first:length(dflp)
        lp(prcId) = plot(dflp(prcId), quantity(prcId), 'o', MarkerSize=ms(prcId), LineWidth=2, Color=cl(prcId, :));hold on;
        lpl(prcId) = string(sprintf('#%d', prcId));
    end
    xlim([0 inf]); ylim([0 inf]); legend(lp, lpl, 'Location','best');
    title(sprintf('%s to days since last paracentesis', quantityTitle));
    xlabel('Days from last paracentesis');ylabel(quantityTitle);
end