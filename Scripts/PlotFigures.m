savefigs  = 1;
%% Plots figures for the paper
color_schema;

addpath('c:\Program Files\Dymola 2021\Mfiles\dymtools\')

datafile = 'HVPGShuntsComparison_extended.mat';
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

shunt_d_res = dymget(dl, 'ascites_Shunts.shunt.d')*1000;
shunt_d_nres = dymget(dl, 'ascites_ShuntStiff.shunt.d')*1000;

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

hvpg_tips = dymget(dl, 'ascites_TIPS.HVPG')/mmHg2SI;
hvpg_sres_tips = dymget(dl, 'ascites_Shunts_TIPS_acute.HVPG')/mmHg2SI;
hvpg_snres_tips = dymget(dl, 'ascites_ShuntStiff_TIPS_acute.HVPG')/mmHg2SI;

q_liver_tips = dymget(dl, 'ascites_TIPS.Q_liver')/L_min2SI;
q_liver_sres = dymget(dl, 'ascites_Shunts.Q_liver')/L_min2SI;
q_liver_snres = dymget(dl, 'ascites_ShuntStiff.Q_liver')/L_min2SI;
q_liver_sres_tips = dymget(dl, 'ascites_Shunts_TIPS_acute.Q_liver')/L_min2SI;
q_liver_snres_tips = dymget(dl, 'ascites_ShuntStiff_TIPS_acute.Q_liver')/L_min2SI;

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
title('Shunt flow');
xlabel('PSG (mmHg)');
ylabel('Shunt flow (ml/min)');
ylim([0, 0.5])
legend('Sensitive', 'Insensitive', 'location', 'northwest')

subplot(222);cla;hold on;
plot(hvpg_sres(inds), shunt_d_res(inds), 'o', 'Linewidth', 2, 'MarkerSize', ms);
plot(hvpg_snres(inds), shunt_d_nres(inds), 's', 'Linewidth', 2, 'MarkerSize', ms);
title('Shunt diameter');
xlabel('PSG (mmHg)');
ylabel('Shunt diameter (mm)');

subplot(223);cla;hold on;
title('PSG')
p_hvpg_ns = plot(R_liver(inds), hvpg_ns(inds), 'd', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms);
p_hvpg_res = plot(R_liver(inds), hvpg_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
p_hvpg_nres = plot(R_liver(inds), hvpg_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('Pressure (mmHg)');
ylim([5 35]);
xlim([5 35]);
xlabel('Liver resistance (mmHg.L/min)');
plot([1 1]*30, [5 35], 'Color', [251 183 74]/255, 'linewidth', 3)
text(30.5, 29, num2str(round(hvpg_ns(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)
text(30.5, 22, num2str(round(hvpg_snres(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)
text(30.5, 18, num2str(round(hvpg_sres(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)
legend([p_hvpg_ns, p_hvpg_res, p_hvpg_nres], 'No shunt', 'shunt, Sensitive', 'shunt, Insensitive', 'location', 'northwest')


subplot(224);cla;hold on;
title('PPV')
plot(R_liver(inds), ppv_ns(inds), 'd', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms);
plot(R_liver(inds), ppv_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
plot(R_liver(inds), ppv_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
ylabel('Pressure (mmHg)');
ylim([10 55]);
xlim([5 35]);
xlabel('Liver resistance (mmHg.L/min)');
plot([1 1]*30, [5 60], 'Color', [251 183 74]/255, 'linewidth', 3)
text(30.5, 52, num2str(round(ppv_ns(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)
text(30.5, 35, num2str(round(ppv_snres(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)
text(30.5, 26, num2str(round(ppv_sres(find(R_liver >= 30, 1, 'first')))), 'fontweight','bold','fontsize',14)

if savefigs
    saveas(f2, 'Figure2_3.png');
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
text(16.8, 24, '2', 'fontweight','bold','fontsize',14)
text(19.8, 24, '3', 'fontweight','bold','fontsize',14)
text(24, 24, '4', 'fontweight','bold','fontsize',14)
yyaxis right;
plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
p_q = plot(R_liver(inds), shunt_q_res(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
ylabel('L/min');
g = gca; g.YAxis(2).Color = color_m;
ylim([0, 0.4])
legend([p_hvpg, p_vol, p_q], 'PSG', 'V_A', 'Q_S', 'position', [0.17, 0.83, 0.1, 0.1])
title('Remodelling sensitive')
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
text(17.6, 24, '3', 'fontweight','bold','fontsize',14)
text(24, 24, '4', 'fontweight','bold','fontsize',14)
yyaxis right;
plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
plot(R_liver(inds), shunt_q_nres(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
ylabel('L/min');
g = gca; g.YAxis(2).Color = color_m;
ylim([0, 0.4])
title('Remodelling insensitive');
xlabel('Liver resistance (mmHg.L/min)')

if savefigs
    saveas(f4, 'Figure4.png');
end
%% Cirrhotic stages plot on PSG
% figure(3);clf;hold on;
% set(gcf, 'Position', [440  360  1000  400])
% 
% % sensitive
% p1_sres = phases_sres == 1;
% p2_sres = phases_sres == 2;
% p3_sres = phases_sres == 3;
% p4_sres = phases_sres == 4;
% p1_sres = find(p2_sres, 1, 'first');
% p2_sres = find(p3_sres, 1, 'first');
% p3_sres = find(p4_sres, 1, 'first');
% p4_sres = length(time);
% 
% h = 25;
% 
% subplot(211);cla;hold on;
% rectangle('Position', [hvpg_sres(1), 0, hvpg_sres(p1_sres) - hvpg_sres(1), h], 'FaceColor', [223, 244, 218]./255, 'Edgecolor', 'None')
% rectangle('Position', [hvpg_sres(p1_sres), 0, hvpg_sres(p2_sres) - hvpg_sres(p1_sres), h], 'FaceColor', [207, 221, 251]./255, 'Edgecolor', 'None')
% rectangle('Position', [hvpg_sres(p2_sres), 0, hvpg_sres(p3_sres) - hvpg_sres(p2_sres), h], 'FaceColor', [252, 200, 190]./255, 'Edgecolor', 'None')
% rectangle('Position', [hvpg_sres(p3_sres), 0, hvpg_sres(end) - hvpg_sres(p3_sres), h], 'FaceColor', [255 161 136]/255, 'Edgecolor', 'None')
% % plot(R_liver, phases_sres);
% plot([0, 40], [1 1]*17, '--', 'Color', color_b, 'Linewidth', 0.5);
% plot([0, 40], [1 1]*5, '--', 'Color', color_r, 'Linewidth', 0.5);
% 
% p_hvpg = plot(hvpg_sres(inds), hvpg_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
% p_vol = plot(hvpg_sres(inds), vol_sres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
% xlim([5, 20]);
% ylim([0, 25]);
% ylabel('mmHg, L');
% text(10, 24, '1', 'fontweight','bold','fontsize',14)
% text(16.8, 24, '2', 'fontweight','bold','fontsize',14)
% text(19.8, 24, '3', 'fontweight','bold','fontsize',14)
% text(27, 24, '4', 'fontweight','bold','fontsize',14)
% yyaxis right;
% plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
% p_q = plot(R_liver(inds), shunt_q_res(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
% ylabel('L/min');
% g = gca; g.YAxis(2).Color = color_m;
% ylim([0, 0.4])
% legend([p_hvpg, p_vol, p_q], 'PSG', 'V_A', 'Q_S', 'position', [0.15, 0.7, 0.1, 0.1])
% title('Remodelling sensitive')
% xlabel('PSG (mmHg)');
% 
% % nonsensitive
% p1_snres = phases_snres == 1;
% p2_snres = phases_snres == 2;
% p3_snres = phases_snres == 3;
% p4_snres = phases_snres == 4;
% p1_snres = find(p2_snres, 1, 'first');
% p2_snres = find(p3_snres, 1, 'first');
% p3_snres = find(p4_snres, 1, 'first');
% p4_snres = length(time);
% 
% h = 25;
% 
% subplot(212);cla;hold on;
% rectangle('Position', [hvpg_snres(1), 0, hvpg_snres(p1_snres) - hvpg_snres(1), h], 'FaceColor', [223, 244, 218]./255, 'Edgecolor', 'None')
% % rectangle('Position', [R_liver(p1_snres), 0, R_liver(p2_snres) - R_liver(p1_snres), h], 'FaceColor', [207, 221, 251]./255, 'Edgecolor', 'None')
% rectangle('Position', [hvpg_snres(p2_snres), 0, hvpg_snres(p3_snres) - hvpg_snres(p2_snres), h], 'FaceColor', [252, 200, 190]./255, 'Edgecolor', 'None')
% rectangle('Position', [hvpg_snres(p3_snres), 0, hvpg_snres(end) - hvpg_snres(p3_snres), h], 'FaceColor', [255 161 136]/255, 'Edgecolor', 'None')
% % plot(R_liver, phases_sres);
% plot([0, 40], [1 1]*17, '--', 'Color', color_b, 'Linewidth', 0.5);
% plot([0, 40], [1 1]*5, '--', 'Color', color_r, 'Linewidth', 0.5);
% plot(hvpg_snres(inds), hvpg_snres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms);
% plot(hvpg_snres(inds), vol_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms);
% xlim([5, 20]);
% ylim([0, 25]);
% ylabel('mmHg, L');
% text(11, 24, '1', 'fontweight','bold','fontsize',14)
% text(17.3, 24, '3', 'fontweight','bold','fontsize',14)
% text(25, 24, '4', 'fontweight','bold','fontsize',14)
% yyaxis right;
% plot([0, 40], [1 1]*0.1, '--', 'Color', color_m, 'Linewidth', 0.5);
% plot(hvpg_snres(inds), shunt_q_nres(inds), 'o', 'Color', color_m, 'Linewidth', ms/2, 'MarkerSize', ms/2);
% ylabel('L/min');
% g = gca; g.YAxis(2).Color = color_m;
% ylim([0, 0.4])
% title('Remodelling insensitive');
% xlabel('PSG (mmHg)');

%% TIPS plotting

f5 = figure(5);clf;
set(gcf, 'Position', [440  360  1000  400])
subplot(131);cla;hold on;
plot(R_liver(inds), - hvpg_sres(inds) + hvpg_sres_tips(inds), 'o', 'Color', color_b, 'Linewidth', 1, 'MarkerSize', ms*0.8, 'MarkerFaceColor', color_b)
plot(R_liver(inds), - hvpg_snres(inds) + hvpg_snres_tips(inds), 's', 'Color', color_r, 'Linewidth', 1, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
plot(R_liver(inds), - hvpg_ns(inds) + hvpg_tips(inds), '+', 'Color', color_s, 'Linewidth', ms/2, 'MarkerSize', ms)
xlim([14, 35]);
title('PSG drop after TIPS placement');
legend('PSG sensitive', 'PSG insensitive', 'No shunts', 'location', 'southwest');
ylabel('(mmHg)');
xlabel('Liver resistance (mmHg.L/min)')

%%
subplot(132);cla;hold on;
plot(R_liver(inds),  R_liver(inds).*0 +1, '+', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms)
plot(R_liver(inds), q_liver_sres(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms)
plot(R_liver(inds), q_liver_snres(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms)
plot(R_liver(inds), q_liver_sres_tips(inds), 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
plot(R_liver(inds), q_liver_snres_tips(inds), 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
plot(R_liver(inds), q_liver_tips(inds), '+', 'Color', color_s, 'Linewidth', ms/2, 'MarkerSize', ms, 'MarkerFaceColor', color_s)
xlim([14, 35]);
ylim([0.1, 1.02])
title('Liver flow before and after TIPS');
xlabel('Liver resistance (mmHg.L/min)');
ylabel('(L/min)');
% legend('No shunt, before TIPS', 'Sensitive, Before TIPS', 'Insensitive, Before TIPS', 'No shunt, after TIPS', 'Sensitive, after TIPS','Insensitive, after TIPS')

% fake legend
axes('position',[0.45 0.4 0.15 0.1]); cla;hold on;
d = 1;
plot(d, 1, '+', 'Color', color_s, 'Linewidth', 2, 'MarkerSize', ms)
plot(2*d, 1, 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms)
plot(3*d, 1, 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms)
text(4*d, 1, 'Before TIPS');

plot(1*d, 0, '+', 'Color', color_s, 'Linewidth', ms/2, 'MarkerSize', ms)
plot(2*d, 0, 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_b)
plot(3*d, 0, 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_r)
text(4*d, 0, 'After TIPS');
xlim([0, 10])
ylim([-0.5, 1.5]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
box on;


%
subplot(133);cla;hold on;
plot(R_liver(inds), -100 + q_liver_tips(inds)./1*100, '+', 'Color', color_s, 'Linewidth', ms/2, 'MarkerSize', ms);
plot(R_liver(inds), -100 + q_liver_snres_tips(inds)./q_liver_snres(inds)*100, 's', 'Color', color_r, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_r);
plot(R_liver(inds), -100 + q_liver_sres_tips(inds)./q_liver_sres(inds)*100, 'o', 'Color', color_b, 'Linewidth', 2, 'MarkerSize', ms, 'MarkerFaceColor', color_b);
xlim([14, 35]);
title('Liver flow decrease after TIPS');
xlabel('Liver resistance (mmHg.L/min)');
ylabel('(%)');
% legend('No shunts', 'remodelling insensitive', 'remodelling sensitive')

if savefigs
    saveas(f5, 'Figure5.png');
end