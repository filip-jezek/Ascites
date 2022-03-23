file = "post_tips_us_14d_ascites.csv";

%     r = csvread("../Data/" + file);
d = LoadTipsCsv("../Data/" + file);

%%
figure(1);clf;hold on;title('pre (dashed) and post TIPS (full) pressures');
colororder('default');
C = colororder;
colororder([C;0.4, 0.4, 0.4;0,0,0;1,0,0;0,1,0;0,0,1]);
set(gca,'ColorOrderIndex',1)
for i = 1:size(d.PPV_base, 1)
    plot([0, 1], [d.PPV_base(i), d.PRA_base(i)], 'o--', 'Linewidth', 1)
end
ylabel('Portal vein (mmHg)');

set(gca,'ColorOrderIndex',1)
for i = 1:size(d.PPV_base, 1)
    plot([0, 1], [d.PPV_TIPS(i), d.PRA_TIPS(i)], 'x-', 'Linewidth', 1)
end

yl = ylim;
yyaxis right;
ylim(yl);
ylabel('Right Atrium (mmHg)');
lab_ppv = "PPV: " + round(nanmean(d.PPV_base)) + " to " + round(nanmean(d.PPV_TIPS));
lab_pra = "PRA: " + round(nanmean(d.PRA_base)) + " to " + round(nanmean(d.PRA_TIPS));
set(gca,'xtick',[0,1],'xticklabel',{lab_ppv, lab_pra})

%%
dp = (-d.PRA_TIPS + d.PPV_TIPS)*133.32;
r = d.StentInternalDiameter/2*1e-2;
L = (d.GraftLinedLength + d.GraftUnlinedLength)*1e-2;
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

figure(2);clf;hold on;
% plot(rand(10, 1), q_est_lpm, 'x');
plot(rand(10, 1), Qmax_lpm, 'o');
plot([0, 1], [1, 1], 'Linewidth', 1);

