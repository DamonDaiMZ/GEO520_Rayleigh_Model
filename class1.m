% P is the product. N is the substrate, i denotes the rare isotope
P = zeros(1,1000);
Pi = zeros(1,1000);
dP = zeros(1,1000);
dPi = zeros(1,1000);
N = zeros(1,1000);
Ni = zeros(1,1000);
Pdelta = zeros(1,1000);
dPdelta = zeros(1,1000);
Ndelta = zeros(1,1000);
R_ref = 0.004;
delta_15 = 5;
Ndelta(1) = delta_15;
N(1) = 1e2;
Ni(1) = ( delta_15/1000 + 1 )*R_ref*N(1);

e = 20;
alpha = 1 - e/1000;
k = 1e-2;
ki = k*alpha;

dt = 1;
for I = 2:1000
    P(I) = P(I-1) + k.*N(I-1).*dt;
    dP(I) = k.*N(I-1).*dt;
    N(I) = N(I-1) - ( P(I) - P(I-1) );
    Pi(I) = Pi(I-1) + ki.*Ni(I-1).*dt;
    dPi(I) = ki.*Ni(I-1).*dt;
    Ni(I) = Ni(I-1) - ( Pi(I) - Pi(I-1) );
    Pdelta(I) = ( Pi(I)./P(I)./R_ref - 1).*1e3;
    dPdelta(I) = ( dPi(I)./dP(I)./R_ref - 1).*1e3;
    Ndelta(I) = ( Ni(I)./N(I)./R_ref - 1).*1e3;
end

figure(1);
logN = log(N./N(1));
plot(logN,Ndelta,"LineWidth",2); hold on;
% plot(logN,Pdelta,"LineWidth",2); hold on;
% plot(logN,dPdelta,"LineWidth",2); hold on;
xlimit = get(gca,'Xlim');
ylimit = get(gca,'Ylim');
yl = [ylimit(2) ylimit(1)];
line(xlimit,yl,'Color','black','LineStyle',':','LineWidth',1); hold on;
slope = ( Ndelta(1) - Ndelta(end) )/( logN(1) - logN(end));
txt = ['slope =' num2str(slope)];
t = text(-3,225,txt);
t.FontSize = 14;
xlabel('ln(f)');
ylabel('substrate \delta^{15}N');
set(gca,'FontSize',14,'linewidth',1);
print(gcf,'-r600','-djpeg',['vs_ln(f)','.jpeg']);

figure(2);
plot(N(2:end-1)./N(1),Ndelta(2:end-1),"LineWidth",2); hold on;
plot(N(2:end-1)./N(1),Pdelta(2:end-1),"LineWidth",2); hold on;
plot(N(2:end-1)./N(1),dPdelta(2:end-1),"LineWidth",2); hold on;
xlimit = get(gca,'Xlim');
line(xlimit,[5 5],'Color','black','LineStyle',':','LineWidth',1); hold on;
legend('substrate','integral product','instaneous product','');
xlabel('f');
ylabel('\delta^{15}N');
set(gca,'FontSize',14,'linewidth',1);
print(gcf,'-r600','-djpeg',['vs_f','.jpeg']);


% substrate, instan. product, integr. product, vs f and log(f)
% one loss term, make it as a steady state.(not due by next next Thurs)
    