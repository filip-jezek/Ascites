savefigs  = 1;
%% Plots figures for the paper
color_schema;

addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = 'HVPGShuntsComparison_IncreasingInflow.mat';
dl = dymload(datafile);

%%
time = dymget(dl, 'Time');

HVPG_nom_max = max(dymget(dl, 'ascites_NoShunts.HVPG'))/mmHg2SI;

mmHgMin_L2SI = (1e+3)*(133.322387415)*60;
dyns_cm52SI = 1e5;

%%
inflow = dymget(dl, 'ramp.y')/L_min2SI;

R_liver = dymget(dl, 'ascites_Shunts.Liver.resistance')/mmHgMin_L2SI;

shunt_q_res = dymget(dl, 'ascites_Shunts.Q_shunt')/L_min2SI;
shunt_q_nres = dymget(dl, 'ascites_ShuntStiff.Q_shunt')/L_min2SI;

shunt_d_res = dymget(dl, 'ascites_Shunts.splenorenalShunt.d')*1000;
shunt_d_nres = dymget(dl, 'ascites_ShuntStiff.splenorenalShunt.d')*1000;

hvpg_ns = dymget(dl, 'ascites_NoShunts.HVPG')/mmHg2SI;
hvpg_sres = dymget(dl, 'ascites_Shunts.HVPG')/mmHg2SI;
hvpg_snres= dymget(dl, 'ascites_ShuntStiff.HVPG')/mmHg2SI;

ppv_ns = dymget(dl, 'ascites_NoShunts.PPV')/mmHg2SI;
ppv_sres = dymget(dl, 'ascites_Shunts.PPV')/mmHg2SI;
ppv_snres= dymget(dl, 'ascites_ShuntStiff.PPV')/mmHg2SI;

vol_sres = dymget(dl, 'ascites_Shunts.levittCase1SsSiIo.Av')*1000;
vol_snres = dymget(dl, 'ascites_ShuntStiff.levittCase1SsSiIo.Av')*1000;

phases_sres = dymget(dl, 'ascites_Shunts.phase');
phases_snres= dymget(dl, 'ascites_ShuntStiff.phase');

q_liver_sres = dymget(dl, 'ascites_Shunts.Q_liver')/L_min2SI;
q_liver_snres = dymget(dl, 'ascites_ShuntStiff.Q_liver')/L_min2SI;

%%
% clf;
% plot(diff(dymget(dl, 'ascites_ShuntStiff.phase')))
%% decimate to individual time points
timepoints = 1:2:max(time);
[~, inds] = min(abs(time - timepoints));
times = time(inds);

%% Figure panel A
f2 = figure(2);clf;
set(gcf, 'Position', [440  10  1000  800])
ms = 8;

subplot(221);cla;hold on;
plot(hvpg_sres(inds), shunt_q_res(inds), 'o', 'Linewidth', 2, 'MarkerSize', ms);
plot(hvpg_snres(inds), shunt_q_nres(inds), 's', 'Linewidth', 2, 'MarkerSize', ms);
title('Shunt flow with increasing splanchnic inflow');
xlabel('PSG (mmHg)');
ylabel('Shunt flow (ml/min)');
ylim([0, 1.5])
xlim([5,35])
cm = get(gca,'colororder');
yyaxis right;

plot(hvpg_sres, inflow, '-', 'Linewidth', 2, 'MarkerSize', ms, 'Color', cm(1,:));
plot(hvpg_snres, inflow, '-', 'Linewidth', 2, 'MarkerSize', ms, 'Color', cm(2,:));
ylabel('*L/min')
% legend('Sensitive', 'Insensitive', 'Splanchnic flow, sens*', 'Splanchnic flow, Insens*', 'location', 'northwest')
set(gca, 'yColor', color_m)
ylim([0, 3])


subplot(222);cla;hold on;
plot(hvpg_sres(inds), shunt_d_res(inds), 'o', 'Linewidth', 2, 'MarkerSize', ms);
plot(hvpg_snres(inds), shunt_d_nres(inds), 's', 'Linewidth', 2, 'MarkerSize', ms);
title('Shunt diameter with increasing splanchnic inflow');
xlabel('PSG (mmHg)');
ylabel('Shunt diameter (mm)');
ylim([0,6])
xlim([5,35])
yyaxis right;
plot(hvpg_sres, inflow, '-', 'Linewidth', 2, 'MarkerSize', ms, 'Color', cm(1,:));
plot(hvpg_snres, inflow, '-', 'Linewidth', 2, 'MarkerSize', ms, 'Color', cm(2,:));
ylabel('*L/min')
% legend('Sensitive', 'Insensitive', 'Splanchnic Q, sens*', 'Splanchnic Q, Insens*', 'location', 'northwest')
legend('Sensitive', 'Insensitive', 'Splanchnic flow, sensitive*', 'Splanchnic flow, Insensitive*', 'location', 'southwest')
set(gca, 'yColor', color_m);
ylim([0, 3])

subplot(223);cla;hold on;
title('PSG with increasing splanchnic inflow')
p_hvpg_ns = plot(R_liver(inds), hvpg_ns(inds), '-d', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms);
p_hvpg_res = plot(R_liver(inds), hvpg_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
p_hvpg_nres = plot(R_liver(inds), hvpg_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('Pressure (mmHg)');
ylim([0 35]);
xlim([5 35]);
xlabel('Liver resistance (mmHg.L/min)');
yyaxis right;

p_inf = plot(R_liver, inflow, '-', 'Color', color_m, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('*L/min');
set(gca, 'yColor', color_m);
ylim([0, 3])

yyaxis left
plot([1 1]*R_liver(find(hvpg_ns >= 30, 1, 'first')), [0 35], 'Color', [251 183 74]/255, 'linewidth', 3)
text(18.5, 29, num2str(round(hvpg_ns(find(hvpg_ns >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_s)
text(18.6, 25, num2str(round(hvpg_snres(find(hvpg_ns >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_r)
text(18.5, 18, num2str(round(hvpg_sres(find(hvpg_ns >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_b)
legend([p_hvpg_ns, p_hvpg_res, p_hvpg_nres, p_inf], 'No shunt', 'shunt, Sensitive', 'shunt, Insensitive', 'Splanchnic inflow*', 'location', 'southeast')


subplot(224);cla;hold on;
title('PPV with increasing splanchnic inflow')
plot(R_liver(inds), ppv_ns(inds), 'd', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms);
plot(R_liver(inds), ppv_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
plot(R_liver(inds), ppv_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('Pressure (mmHg)');
ylim([0 80]);
xlim([5 35]);
xlabel('Liver resistance (mmHg.L/min)');
yyaxis right;

plot(R_liver, inflow, '-', 'Color', color_m, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('*L/min')
set(gca, 'yColor', color_m);
ylim([0, 3])

yyaxis left
plot([1 1]*R_liver(find(ppv_ns >= 55, 1, 'first')), [0 100], 'Color', [251 183 74]/255, 'linewidth', 3)
text(15, 55, num2str(round(ppv_ns(find(ppv_ns >= 55, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_s)
text(15.5, 34.5, num2str(round(ppv_snres(find(ppv_ns >= 55, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_r)
text(17.9, 27, num2str(round(ppv_sres(find(ppv_ns >= 55, 1, 'first')))), 'fontweight','bold','fontsize',14, 'Color', color_b)
legend([p_hvpg_ns, p_hvpg_res, p_hvpg_nres, p_inf], 'No shunt', 'shunt, Sensitive', 'shunt, Insensitive', 'Splanchnic inflow*', 'location', 'southeast')

if savefigs
    saveas(f2, 'Suppl_Figure2_3.png');
end
%% cirrhotic stages plots on liver resistance

f4 = figure(4);clf;hold on;
set(gcf, 'Position', [440  60  500  600])

% sensitive
p1_sres = phases_sres == 1;
p2_sres = phases_sres == 2;
p3_sres = phases_sres == 3;
p4_sres = phases_sres == 4;
p1_sres = find(p2_sres, 1, 'first');
p2_sres = find(p3_sres, 1, 'first');
p3_sres = find(p4_sres, 1, 'first');
p4_sres = length(time);

h = 25;

subplot(211);cla;hold on;
rectangle('Position', [R_liver(1), 0, R_liver(p1_sres) - R_liver(1), h], 'FaceColor', [223, 244, 218]./255, 'Edgecolor', 'None')
rectangle('Position', [R_liver(p1_sres), 0, R_liver(p2_sres) - R_liver(p1_sres), h], 'FaceColor', [207, 221, 251]./255, 'Edgecolor', 'None')
rectangle('Position', [R_liver(p2_sres), 0, R_liver(p3_sres) - R_liver(p2_sres), h], 'FaceColor', [252, 200, 190]./255, 'Edgecolor', 'None')
rectangle('Position', [R_liver(p3_sres), 0, R_liver(end) - R_liver(p3_sres), h], 'FaceColor', [255 161 136]/255, 'Edgecolor', 'None')
% plot(R_liver, phases_sres);
plot([0, 40], [1 1]*17, '--', 'Color', color_b, 'Linewidth', 0.5);
plot([0, 40], [1 1]*5, '--', 'Color', color_r, 'Linewidth', 0.5);

p_hvpg = plot(R_liver(inds), hvpg_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
p_vol = plot(R_liver(inds), vol_sres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
xlim([5, 25]);
ylim([0, 25]);
ylabel('mmHg, L');
text(11, 24, '1', 'fontweight','bold','fontsize',14)
text(13.5, 24, '2', 'fontweight','bold','fontsize',14)
text(14.4, 24, '3', 'fontweight','bold','fontsize',14)
text(19.5, 24, '4', 'fontweight','bold','fontsize',14)
yyaxis right;
plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
p_q = plot(R_liver(inds), shunt_q_res(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
ylabel('L/min');
g = gca; g.YAxis(2).Color = color_m;
ylim([0, 0.4])
legend([p_hvpg, p_vol, p_q], 'PSG', 'V_A', 'Q_S', 'position', [0.17, 0.83, 0.1, 0.1])
title('Remodelling sensitive with increasing splanchnic inflow')
% xlabel('Liver resistance (mmHg.L/min)');

% nonsensitive
p1_snres = phases_snres == 1;
p2_snres = phases_snres == 2;
p3_snres = phases_snres == 3;
p4_snres = phases_snres == 4;
p1_snres = find(p1_snres, 1, 'last');
p2_snres = find(p3_snres, 1, 'first');
p3_snres = find(p4_snres, 1, 'first');
p4_snres = length(time);

h = 25;

subplot(212);cla;hold on;
rectangle('Position', [R_liver(1), 0, R_liver(p1_snres) - R_liver(1), h], 'FaceColor', [223, 244, 218]./255, 'Edgecolor', 'None')
% rectangle('Position', [R_liver(p1_snres), 0, R_liver(p2_snres) - R_liver(p1_snres), h], 'FaceColor', [207, 221, 251]./255, 'Edgecolor', 'None')
rectangle('Position', [R_liver(p2_snres), 0, R_liver(p3_snres) - R_liver(p2_snres), h], 'FaceColor', [252, 200, 190]./255, 'Edgecolor', 'None')
rectangle('Position', [R_liver(p3_snres), 0, R_liver(end) - R_liver(p3_snres), h], 'FaceColor', [255 161 136]/255, 'Edgecolor', 'None')
% plot(R_liver, phases_sres);
plot([0, 40], [1 1]*17, '--', 'Color', color_b, 'Linewidth', 0.5);
plot([0, 40], [1 1]*5, '--', 'Color', color_r, 'Linewidth', 0.5);
plot(R_liver(inds), hvpg_snres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
plot(R_liver(inds), vol_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
xlim([5, 25]);
ylim([0, 25]);
ylabel('mmHg, L');
text(11, 24, '1', 'fontweight','bold','fontsize',14)
text(15.6, 24, '3', 'fontweight','bold','fontsize',14)
text(22, 24, '4', 'fontweight','bold','fontsize',14)
yyaxis right;
plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
plot(R_liver(inds), shunt_q_nres(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
ylabel('L/min');
g = gca; g.YAxis(2).Color = color_m;
ylim([0, 0.4])
title('Remodelling insensitive with increasing splanchnic inflow');
xlabel('Liver resistance (mmHg.L/min)')

if savefigs
    saveas(f4, 'Suppl_Figure4.png');
end
