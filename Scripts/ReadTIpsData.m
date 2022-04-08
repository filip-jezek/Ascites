file = "post_tips_us_14d_evb.csv";
d_evb = LoadTipsCsv("../Data/" + file);

file = "post_tips_us_14d_ascites.csv";
d_asc = LoadTipsCsv("../Data/" + file);

%% ascites figure
figure(1);clf;hold on;title('pre (dashed) and post TIPS (full) pressures EVB');
colororder('default');
C = colororder;
colororder([C;0.4, 0.4, 0.4;0,0,0;1,0,0;0,1,0;0,0,1]);
set(gca,'ColorOrderIndex',1)
for i = 1:size(d_evb.PPV_base, 1)
    plot([0, 1], [d_evb.PPV_base(i), d_evb.PRA_base(i)], 'd-', 'Linewidth', 0.5)
end

ylabel('Portal vein (mmHg)');


set(gca,'ColorOrderIndex',1)
for i = 1:size(d_evb.PPV_base, 1)
    plot([0, 1], [d_evb.PPV_TIPS(i), d_evb.PRA_TIPS(i)], 's-', 'Linewidth', 2)
end


yl = [0, 40];
ylim(yl);
yyaxis right;
ylim(yl);
ylabel('Right Atrium (mmHg)');
lab_ppv = "PPV: " + round(nanmean(d_evb.PPV_base)) + " to " + round(nanmean(d_evb.PPV_TIPS));
lab_pra = "PRA: " + round(nanmean(d_evb.PRA_base)) + " to " + round(nanmean(d_evb.PRA_TIPS));
set(gca,'xtick',[0,1],'xticklabel',{lab_ppv, lab_pra})

%% EVB figure
figure(2);clf;hold on;title('pre (dashed) and post TIPS (full) pressures ASCITES');
colororder('default');
C = colororder;
colororder([C;0.4, 0.4, 0.4;0,0,0;1,0,0;0,1,0;0,0,1]);
set(gca,'ColorOrderIndex',1)
for i = 1:size(d_asc.PPV_base, 1)
    plot([0, 1], [d_asc.PPV_base(i), d_asc.PRA_base(i)], 'o-', 'Linewidth', 0.5)
end

ylabel('Portal vein (mmHg)');


set(gca,'ColorOrderIndex',1)
for i = 1:size(d_asc.PPV_base, 1)
    plot([0, 1], [d_asc.PPV_TIPS(i), d_asc.PRA_TIPS(i)], 'x-', 'Linewidth', 2)
end


yl = [0, 40];
ylim(yl);
yyaxis right;
ylim(yl);
ylabel('Right Atrium (mmHg)');
lab_ppv = "PPV: " + round(nanmean(d_asc.PPV_base)) + " to " + round(nanmean(d_asc.PPV_TIPS));
lab_pra = "PRA: " + round(nanmean(d_asc.PRA_base)) + " to " + round(nanmean(d_asc.PRA_TIPS));
set(gca,'xtick',[0,1],'xticklabel',{lab_ppv, lab_pra})
%%
dp = (-d_asc.PRA_TIPS + d_asc.PPV_TIPS)*133.32;
r = d_asc.StentInternalDiameter/2*1e-2;
L = (d_asc.GraftLinedLength + d_asc.GraftUnlinedLength)*1e-2;
u = 6e-3;
ro = 1e3;
q_estSI = pi*dp.*(r.^4)./(8*L*u);
q_est_lpm = q_estSI*6e4;
%Check with Bernou
QmaxSI = pi*r.^2.*sqrt(2*dp/ro);
Qmax_lpm = QmaxSI*6e4;

% check Reynolds
Q = 1/1e4;
Re = ro*Q*2*r./(u*pi*r.^2);
L./r
Re/48

figure(3);clf;hold on;
% plot(rand(10, 1), q_est_lpm, 'x');
plot(rand(10, 1), Qmax_lpm, 'o');
plot([0, 1], [1, 1], 'Linewidth', 1);

