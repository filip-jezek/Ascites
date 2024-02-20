% Load paracentesis
dataTable = readtable('../data/pressure-volume-time_Adj.csv')

plot(dataTable.Var1, dataTable.Liters_Drained__, dataTable.Var1,dataTable.Record_ID_Adj);
%%
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
end
%% squish empty procedure    
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
    end
end
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

    nexttile;
    
    for prcId = 1:length(patients{patId}.prc)
        prc = patients{patId}.prc{prcId};
        plot(prc.Drained, prc.Pressure, 'o-');hold on;
        
    end
    title(sprintf('Patient %d', patId));
    xlim([0 inf]);ylim([0 inf]);
end
%% single procedure MODEL fit and estimation
figure(2);clf;
tiledlayout('flow');
nexttile;

% linear model
pm = @(a, b, c, x) max(0, a*(x -b) + c);
pm_0 = @(a, b, c, x) pm(a, 0, c, x) + 0*b;

patId = 8;
% patId = 2;
patId = 11;
% prcId = 1;
x_ = 0:10;
hold on;
cl = lines(length(patients{patId}.prc));
for prcId = 1:size(patients{patId}.prc, 2)
    if prcId == 1
        lw = 3;
    else
        lw = 1;
    end

    [ae gd] = fit(patients{patId}.prc{prcId}.Drained', patients{patId}.prc{prcId}.Pressure', pm_0, ...
        'StartPoint', [-1, 0, patients{patId}.prc{prcId}.Pressure(1)]);
    plot(patients{patId}.prc{prcId}.Drained, patients{patId}.prc{prcId}.Pressure, 'o-', 'Color', cl(prcId, :), LineWidth=lw);
    plot(x_, pm(ae.a, ae.b, ae.c, x_), '--', LineWidth=lw, Color=cl(prcId, :));
end
title(sprintf('Single procedure best fit: %0.2fx + %0.2f', ae.a, ae.c))
nexttile;
% Fit all procedures in single patient
fitPat = @(params) evalFit(params, patients{patId}.prc, false)
fr = fminsearch(fitPat, [-1 10], options);
%
% cla;hold on;
[costAvg b costs] = evalFit([fr], patients{patId}.prc, true);
title(sprintf('All procedures best fit: %0.2f(x - V_0) + %0.2f', fr(1), fr(2)));
% Visualize date shift
% max drained volume, which is empty ambdomen
V_dmax = fsolve(@(x)pm(fr(1), 0, fr(2), x), 1e-6);
patients{patId}.V_dmax = V_dmax;

for prcId = 1:length(patients{patId}.prc)
    patients{patId}.dayOfProc(prcId) = patients{patId}.prc{prcId}.Date;
    patients{patId}.b(prcId) = b(prcId); % shift
    % estimate opening volume by openinig pressure, using identified
    % pressure-volume model
    p_open = patients{patId}.prc{prcId}.Pressure(1);
    v_open = fsolve(@(x)pm(fr(1), 0, fr(2), x) - p_open, 1e-6);
    patients{patId}.v_open(prcId) = V_dmax - v_open;
    plot([patients{patId}.prc{prcId}.Drained(1) - b(prcId), v_open], [p_open p_open], 's--', Color=cl(prcId, :), MarkerFaceColor=cl(prcId, :))
end
%
nexttile;hold on;
% clf;hold on;
dflp_inf = 360; % days from last paracentesis - considering steady state
% days from last paracentesis
dflp = days(diff(patients{patId}.dayOfProc))
d_dflp = [dflp_inf dflp]; % days from last paracentesis DIFF
% max out at dflp_inf
d_dflp(d_dflp > dflp_inf) = dflp_inf;

vols = patients{patId}.v_open(1:end);
plot(d_dflp(2:end), vols(2:end), 'k--');
% first paracentesis
plot([0 max(d_dflp(2:end))*1.5], [vols(1) vols(1)], '--', LineWidth=2, Color=cl(1, :))
for prcId = 2:length(patients{patId}.dayOfProc)
    plot(d_dflp(prcId), vols(prcId), 'o', MarkerSize=1/costs(prcId), LineWidth=2, Color=cl(prcId, :))
end
xlim([0 inf]);
ylim([0 inf]);

% Ascites Generation Curve
% agc = a*(1- exp(t*b))

%% Identify all and plot
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.01, 'MaxIter', 50);

figure(3);clf;
tiledlayout('flow');
params = zeros(length(patId), 2);
for patId = 1:min(120,length(patients))
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end

    nexttile;
    prc = patients{patId}.prc;
    fitPat = @(params) evalFit(params, prc, false);
    fr = fminsearch(fitPat, [-1 10], options);
    [costAvg, shift, costs] = evalFit(fr, prc, true);
    patients{patId}.params = fr;
    patients{patId}.shifts = shift;
    patients{patId}.costs = costs;

    % max drained volume, which is zero-pressure volume
    % we get it as a drained volume, at which the model crosses 0 pressure
    % (from above, thus smallest positive number will do)
    V_dmax = fsolve(@(x)pm(fr(1), 0, fr(2), x), 1e-6);
    patients{patId}.V_dmax = V_dmax;

    for prcId = 1:length(patients{patId}.prc)
        p_open = patients{patId}.prc{prcId}.Pressure(1);
        v_at_open = fsolve(@(x)pm(fr(1), 0, fr(2), x) - p_open, 1e-6);
        patients{patId}.p_open(prcId) = p_open;
        patients{patId}.v_open(prcId) = V_dmax - v_at_open;
    end

    % params(patId, :) = fr;
    title(sprintf('Patient %d costs %0.2f', patId, costAvg));
    xlim([-inf inf]);ylim([0 inf]);
end

%% estimate the total volume
figure(5);clf;tl = tiledlayout('flow');
title(tl,'Opening volume by days of last paracentesis')
% hold on;
% colororder(jet(28));
% colororder(gca, "reef")
for patId = 1:min(120, length(patients))
    if isempty(patients{patId})
        continue;
    elseif length(patients{patId}.prc) == 1
        % ignore single procedure patients
        continue;
    end
    nexttile; hold on;
    cl = lines(length(patients{patId}.prc));
    dflp_inf = 360; % days from last paracentesis - considering steady state
    % days from last paracentesis
    dflp = days(diff(patients{patId}.dayOfPrc));
    d_dflp = [dflp_inf dflp]; % days from last paracentesis DIFF
    % max out at dflp_inf
    d_dflp(d_dflp > dflp_inf) = dflp_inf;

    vols = patients{patId}.v_open;
    plot(d_dflp(2:end), vols(2:end), '--');
    % first paracentesis
    plot([0 max(d_dflp(2:end))*1.5], [vols(1) vols(1)], '--', LineWidth=2, Color=cl(1, :))
    for prcId = 2:length(patients{patId}.dayOfPrc)
        ms = 1/patients{patId}.costs(prcId);
        if isnan(ms)
            % single measurement - not enough params to fit curve
            plot(d_dflp(prcId), vols(prcId), 'x', MarkerSize=5, LineWidth=2, Color=cl(prcId, :))
        else
            plot(d_dflp(prcId), vols(prcId), 'o', MarkerSize=ms, LineWidth=2, Color=cl(prcId, :))
        end        
    end
    xlim([0 60]);
    ylim([0 30]);
    title(sprintf('Patient %d', patId));    
end
%%
figure(6);clf;
tl = tiledlayout('flow');
title(tl, 'Some correlations')

ax1 = nexttile;hold on;
xlabel('slope (a)');ylabel('offset (c)');

ax2 = nexttile;hold on;
xlabel('Opening volume');ylabel('slope (a)');

ax3 = nexttile;hold on;
xlabel('Opening volume');ylabel('offset (c)');

ax4 = nexttile;hold on;
xlabel('Opening volume');ylabel('shift (b)');

ax5 = nexttile;hold on;
xlabel('Opening volume');ylabel('Openinig pressure');


for patId = 1:length(patients)
    if ~isfield(patients{patId}, 'params')
        continue;
    end
    plot(ax1, patients{patId}.params(1), patients{patId}.params(2), 'o')
    plot(ax2, patients{patId}.v_open, patients{patId}.params(1), 'o')
    plot(ax3, patients{patId}.v_open, patients{patId}.params(2), 'o')
    plot(ax4, patients{patId}.v_open, patients{patId}.shifts, 'o')
    plot(ax5, patients{patId}.v_open, patients{patId}.p_open, 'o')
end
ylim([0 inf]);xlim([0 inf])
%%
function [costAvg shift costs]= evalFit(params, prcSet, showPlots)
    useReference = false;
    a = params(1);
    c = params(2);
    pm_base = @(b, x) max(0, a*(x - b) + c);
    pm_0 = @(b, x) pm_base(0, x) + 0*b;
    costs = zeros(length(prcSet), 1);
    shift = [];
    
    % reference fit - just for the costs here
    % [ae gd] = fit(prcSet{1}.Drained', prcSet{1}.Pressure', pm, 'StartPoint', [0]);
    % cost = gd.rmse;
    
    if showPlots
    %     plot(prcSet{1}.Drained, prcSet{1}.Pressure);
    % 
    %     hold on;
        x_ = 0:10;
    %     plot(x_, pm_base(0, x_), '--')
    end
    cl = lines(length(prcSet));
    for i = 1:length(prcSet)
        if useReference && i == 1
            % all shifts relative to first one
            pm = pm_0;
        else
            pm = pm_base;
        end
        % shift other drains
        [ae gd] = fit(prcSet{i}.Drained', prcSet{i}.Pressure', pm, 'StartPoint', [0]);        
        costs(i) = gd.rmse/(length(prcSet{i}.Drained));
        shift(i) = ae.b;
        if showPlots
            plot(prcSet{i}.Drained-ae.b, prcSet{i}.Pressure, 'o-', Color=cl(i, :));        
            hold on;
            plot(x_ - ae.b, pm(ae.b, x_), '--', LineWidth=2, Color=cl(i, :));
            leg{i} = sprintf('V_{0} %0.2f', shift(i));
        end
    end
    costAvg = mean(costs);
    if showPlots
        legend(leg)
    end
end


