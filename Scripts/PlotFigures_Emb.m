% plots the prediction of the splenic a. embolization on ascites and PSG

savefigs  = 0;
%% Plots figures for the paper
color_schema;

addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = 'EMB-HVPGShuntsComparison.mat';
dl = dymload(datafile);

%
time = dymget(dl, 'Time');

HVPG_nom_max = max(dymget(dl, 'ascites_NoShunts.HVPG'))/mmHg2SI;

mmHgMin_L2SI = (1e+3)*(133.322387415)*60;
dyns_cm52SI = 1e5;

%
R_liver = dymget(dl, 'ascites_Shunts.Liver.resistance')/mmHgMin_L2SI;

shunt_q_res = dymget(dl, 'ascites_Shunts.Q_shunt')/L_min2SI;
shunt_q_nres = dymget(dl, 'ascites_ShuntStiff.Q_shunt')/L_min2SI;

hvpg_ns = dymget(dl, 'ascites_NoShunts.HVPG')/mmHg2SI;
hvpg_nsEmb = dymget(dl, 'ascites_NoShuntsEmb.HVPG')/mmHg2SI;
hvpg_s = dymget(dl, 'ascites_Shunts.HVPG')/mmHg2SI;
hvpg_sEmb = dymget(dl, 'ascites_ShuntsEmb.HVPG')/mmHg2SI;
%

vol_ns = dymget(dl, 'ascites_NoShunts.levittCase1SsSiIo.Av')*1000;
vol_s = dymget(dl, 'ascites_Shunts.levittCase1SsSiIo.Av')*1000;
vol_nsEmb = dymget(dl, 'ascites_NoShuntsEmb.levittCase1SsSiIo.Av')*1000;
vol_sEmb = dymget(dl, 'ascites_ShuntsEmb.levittCase1SsSiIo.Av')*1000;

%
% clf;
% plot(diff(dymget(dl, 'ascites_ShuntStiff.phase')))
% decimate to individual time points
timepoints = 1:2:max(time);
[~, inds] = min(abs(time - timepoints));
times = time(inds);
ms = 8;
%% EMB plotting - to liver resistance
cursorPos = 25; % at resistance 20
f1 = figure(1);clf;
set(gcf, 'Position', [440  360  1000  400])
subplot(131);cla;hold on;
nsf = plot(R_liver(inds), hvpg_ns(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
plot(R_liver(inds), hvpg_nsEmb(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
sf = plot(R_liver(inds), hvpg_s(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms)
plot(R_liver(inds), hvpg_sEmb(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms,'MarkerFaceColor', color_r)
plot([cursorPos cursorPos], [0 50], '--', 'Color', [251 183 74]/255, 'Linewidth', 1)
xlim([10, 30]);ylim([5, 30]);
title('PSG after SA embolization');
legend('No shunts', 'No shunts emb', 'Shunt', 'Shunt emb', 'location', 'southeast');
ylabel('(mmHg)');
xlabel('Liver resistance (mmHg.L/min)')

%
subplot(132);cla;hold on;
plot(R_liver(inds), vol_ns(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
plot(R_liver(inds), vol_nsEmb(inds), 's', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
plot(R_liver(inds), vol_s(inds), 'o', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms*0.8)
plot(R_liver(inds), vol_sEmb(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
xlim([10, 30]); ylim([0 20]);
plot([cursorPos cursorPos], [0 50], '--', 'Color', [251 183 74]/255, 'Linewidth', 1)
title('Ascites volume after SA embolization');
% legend('No shunts', 'No shunts emb', 'Shunt', 'Shunt emb', 'location', 'southwest');
ylabel('(L)');
xlabel('Liver resistance (mmHg.L/min)')
% legend('No shunt, before TIPS', 'Sensitive, Before TIPS', 'Insensitive, Before TIPS', 'No shunt, after TIPS', 'Sensitive, after TIPS','Insensitive, after TIPS')

%
subplot(133);cla;hold on;
plot(R_liver(inds), vol_nsEmb(inds) - vol_ns(inds), 'o-', 'Color', color_b , 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(R_liver(inds), vol_sEmb(inds) - vol_s(inds), 's-', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
plot([cursorPos cursorPos], [-50 50], '--', 'Color', color_y, 'Linewidth', 1)
xlim([10, 30]);
ylim([-10, 1]);
title('Ascites volume decrease after EMB');
legend('No shunt', 'Shunt', 'Location', 'southwest');
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


%% differences on liver resistance
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

%% EMB plotting to pre x post HVPG
f2 = figure(2);clf;
PPG_change = [21.5, 17.6]% pre-emb, post emb from 10.1002/lt.20762

% Data from doi 10.1002/lt.20762
data = [-0.0012366671819445951, 27.015416180665838;
1.0212341376354368, 20.945347742379745;
0.006661139139110928, 25.97015135119943;
1.021487092286289, 21.552438904425305;
-0.0016582582666984091, 26.00359757725657;
1.0212341376354368, 20.945347742379745;
-0.00997765567250819, 26.037043803313708;
1.0196039854410552, 17.032982475863914;
-0.0033165165333968183, 22.02377773718012;
1.0208687586953167, 20.068438286091713;
-0.010441405865737607, 24.924043339563518;
1.00584606304192, 24.013968717941516;
0.004609395859975551, 21.04596748127433;
1.019449402043312, 16.66198232128052;
0.00456723675150017, 20.944785620933406;
1.0095982236962298, 13.019154288283989;
-0.004173751739063203, 19.96641324358128;
1.0192245534647764, 16.122345732795573;
-0.004581289787658482, 18.98832192695232;
1.0025857586531568, 16.189238184909854;
-0.012928793265785332, 18.954313579448844;
0.9940415126688122, 15.68304782248205;
-0.02254106999817318, 15.884849421717567;
1.0187326971992305, 14.941890695484762;
0.0020939023876107576, 15.00878314759904;
1.0179176211020395, 12.985708062226848];

psg_pre = data(1:2:end, 2);
psg_post = data(2:2:end, 2);

set(gcf, 'Position', [440  360  1000  400])
subplot(131);cla;hold on;
nsf = plot(hvpg_ns(inds), hvpg_nsEmb(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
sf = plot(hvpg_s(inds), hvpg_sEmb(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms)
dataf = plot(psg_pre, psg_post, 'x', 'Color', color_y, 'Linewidth', 2, 'MarkerSize', ms)
xlim([10, 30]);
title('PSG after SA embolization');
legend('No shunts', 'Shunt', 'Data Luca 2006', 'location', 'northwest');
ylabel('post-embolization PSG (mmHg)');
xlabel('Pre-embolization PSG (mmHg)')

%
subplot(132);cla;hold on;
plot(vol_ns(inds), vol_nsEmb(inds),'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8)
% plot(hvpg_ns(inds), vol_nsEmb(inds), 's', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
plot(vol_s(inds), vol_sEmb(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms)
% plot(hvpg_s(inds), vol_sEmb(inds), 'x', 'Color', color_r, 'Linewidth', ms/2, 'MarkerSize', ms)
xlim([0, 30]);
ylim([0 15])
title('Ascites volume after SA embolization');
% legend('No shunts', 'No shunts emb', 'Shunt', 'Shunt emb', 'location', 'southwest');
ylabel('(L)');
xlabel('Pre-embolization ascites volume (mmHg)')
% legend('No shunt, before TIPS', 'Sensitive, Before TIPS', 'Insensitive, Before TIPS', 'No shunt, after TIPS', 'Sensitive, after TIPS','Insensitive, after TIPS')

%
subplot(133);cla;hold on;
plot(hvpg_ns(inds), vol_nsEmb(inds) - vol_ns(inds), 'o-', 'Color', color_b , 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(hvpg_s(inds), vol_sEmb(inds) - vol_s(inds), 's-', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
xlim([10, 30]);
title('Ascites volume decrease after EMB');
legend('No shunt', 'Shunt');
xlabel('Pre-embolization PSG (mmHg)')
ylabel('(L)');

%%
if savefigs
    saveas(f1, 'EMB_Figure1Plus.png');
    saveas(f2, 'EMB_Figure2Diffs.png');
end
