savefigs  = 0;
%% Plots figures for the paper
color_schema;

addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = 'EMB-HVPGShuntsComparison.mat';
dl = dymload(datafile);

%%
time = dymget(dl, 'Time');

HVPG_nom_max = max(dymget(dl, 'ascites_NoShunts.HVPG'))/mmHg2SI;

mmHgMin_L2SI = (1e+3)*(133.322387415)*60;
dyns_cm52SI = 1e5;

%%
R_liver = dymget(dl, 'ascites_Shunts.Liver.resistance')/mmHgMin_L2SI;

shunt_q_res = dymget(dl, 'ascites_Shunts.Q_shunt')/L_min2SI;
shunt_q_nres = dymget(dl, 'ascites_ShuntStiff.Q_shunt')/L_min2SI;

hvpg_ns = dymget(dl, 'ascites_NoShunts.HVPG')/mmHg2SI;
hvpg_nsEmb = dymget(dl, 'ascites_NoShuntsEmb.HVPG')/mmHg2SI;
hvpg_s = dymget(dl, 'ascites_Shunts.HVPG')/mmHg2SI;
hvpg_sEmb = dymget(dl, 'ascites_ShuntsEmb.HVPG')/mmHg2SI;
%%

vol_ns = dymget(dl, 'ascites_NoShunts.levittCase1SsSiIo.Av')*1000;
vol_s = dymget(dl, 'ascites_Shunts.levittCase1SsSiIo.Av')*1000;
vol_nsEmb = dymget(dl, 'ascites_NoShuntsEmb.levittCase1SsSiIo.Av')*1000;
vol_sEmb = dymget(dl, 'ascites_ShuntsEmb.levittCase1SsSiIo.Av')*1000;

%%
% clf;
% plot(diff(dymget(dl, 'ascites_ShuntStiff.phase')))
%% decimate to individual time points
timepoints = 1:2:max(time);
[~, inds] = min(abs(time - timepoints));
times = time(inds);

%% EMB plotting
ms = 8;
f1 = figure(1);clf;
set(gcf, 'Position', [440  360  1000  400])
subplot(131);cla;hold on;
ns = plot(R_liver(inds), hvpg_ns(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
plot(R_liver(inds), hvpg_nsEmb(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
s = plot(R_liver(inds), hvpg_s(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms)
plot(R_liver(inds), hvpg_sEmb(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms,'MarkerFaceColor', color_r)
xlim([10, 35]);
title('PSG after SA embolization');
legend('No shunts', 'No shunts emb', 'Shunt', 'Shunt emb', 'location', 'northwest');
ylabel('(mmHg)');
xlabel('Liver resistance (mmHg.L/min)')

%
subplot(132);cla;hold on;
plot(R_liver(inds), vol_ns(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
plot(R_liver(inds), vol_nsEmb(inds), 's', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
plot(R_liver(inds), vol_s(inds), '+', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms)
plot(R_liver(inds), vol_sEmb(inds), 'x', 'Color', color_r, 'Linewidth', ms/2, 'MarkerSize', ms)
xlim([10, 35]);
ylim([0 25])
title('Ascites volume after SA embolization');
% legend('No shunts', 'No shunts emb', 'Shunt', 'Shunt emb', 'location', 'southwest');
ylabel('(L)');
xlabel('Liver resistance (mmHg.L/min)')
% legend('No shunt, before TIPS', 'Sensitive, Before TIPS', 'Insensitive, Before TIPS', 'No shunt, after TIPS', 'Sensitive, after TIPS','Insensitive, after TIPS')

%
subplot(133);cla;hold on;
plot(R_liver(inds), vol_nsEmb(inds) - vol_ns(inds), 'o-', 'Color', color_b , 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(R_liver(inds), vol_sEmb(inds) - vol_s(inds), 's-', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
xlim([14, 35]);
title('Ascites volume decrease after EMB');
legend('No shunt', 'Shunt');
xlabel('Liver resistance (mmHg.L/min)');
ylabel('(L)');

%%
% fake legend
% axes('position',[0.45 0.4 0.15 0.1]); cla;hold on;
% d = 1;
% plot(d, 1, 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms)
% plot(2*d, 1, 's', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms)
% text(3*d, 1, 'Before emb');
% 
% plot(1*d, 0, 'o', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
% plot(2*d, 0, 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
% text(3*d, 0, 'After emb');
% xlim([0, 10])
% ylim([-0.5, 1.5]);
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% box on;


%%
f2 = figure(2);clf;
subplot(121);cla;hold on;
plot(R_liver(inds), hvpg_nsEmb(inds) - hvpg_ns(inds), 'o-', 'Color', color_b , 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(R_liver(inds), hvpg_sEmb(inds) - hvpg_s(inds), 's-', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)

xlim([4, 35]);
title('HVPG decrease after EMB');
legend('No shunt', 'Shunt');
xlabel('Liver resistance (mmHg.L/min)');
ylabel('(mmHg)');
% legend('No shunts', 'remodelling insensitive', 'remodelling sensitive')

subplot(122);cla;hold on;
plot(R_liver(inds), vol_nsEmb(inds) - vol_ns(inds), 'o-', 'Color', color_b , 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(R_liver(inds), vol_sEmb(inds) - vol_s(inds), 's-', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
xlim([14, 35]);
title('Ascites volume decrease after EMB');
legend('No shunt', 'Shunt');
xlabel('Liver resistance (mmHg.L/min)');
ylabel('(L)');


if savefigs
    saveas(f1, 'EMB_Figure1Plus.png');
    saveas(f2, 'EMB_Figure2Diffs.png');
end