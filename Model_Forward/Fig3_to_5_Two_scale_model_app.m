clc;  clear all;
%%%%%%%%%%%%%%%%%%%%%%%     Sensitivity simulation of rock physics models  %%%%%%%%%%%%%%%%%%%%
Fontsize1 = 14;
load ./model_data./Kmd_gas.mat;
load ./model_data./umd_gas.mat;
Kmd_gas = Kmd;
umd_gas = umd;
load ./model_data./kmd_water.mat;
load ./model_data./umd_water.mat;
Kmd_water = Kmd;
umd_water = umd;
load ./model_data./kmd_oil.mat;
load ./model_data./umd_oil.mat;
Kmd_oil = Kmd;
umd_oil = umd;

load ./model_data./K_dry.mat;
load ./model_data./u_dry.mat;

FrequencyRange=[10.^(-4:0.01:-0.01) 1:1:200 10.^(2.4:0.01:8)];
w=2*pi*FrequencyRange;

Ks = 39e9;
us = 35e9;
rho_s = 2650;
por = 0.12;   
[n1,n2] = size(Kmd_oil);
%%%%%%%%%%%%%%% Parameters of fluids ----------------------------------------------------
Kf_oil     = 4.66*10^(9);    % Glycerol
rho_oil = 1261;            % (Kg/m^3)
vis_oil    = 1400*10^(-3);      % (pa.s) 1pa.s=1000mpa.s=1000cP

Kf_gas  = 0.01*10^(9);              % gas 
rho_gas = 50;
vis_gas = 1*10^(-5);          

Kf_water     = 2.25*10^(9);        % Water
rho_water = 1000;                
vis_water    = 1*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;                         %  Choose 1- water + gas, 2- water + oil, 3- Oil + gas
S_1 = [0.1 0.3 0.5 0.7 0.9];   %  The saturation of fluid 1
perm0 = [1e-16 1e-15 1e-14 1e-13 1e-12] ;

if n ==1
    Kf = [Kf_water Kf_gas];                 
    rho_f = [rho_water rho_gas];            
    vis_f = [vis_water vis_gas];            
    Kmd1 = Kmd_water;   Kmd2 = Kmd_gas;
    umd1 = umd_water;   umd2 = umd_gas;
elseif n == 2
    Kf = [Kf_water Kf_oil];                 
    rho_f = [rho_water rho_oil];            
    vis_f = [vis_water vis_oil];            
    Kmd1 = Kmd_water;   Kmd2 = Kmd_oil;
    umd1 = umd_water;   umd2 = umd_oil;
else
    Kf = [Kf_oil Kf_gas];                
    rho_f = [rho_oil rho_gas];           
    vis_f = [vis_oil vis_gas];            
    Kmd1 = Kmd_oil;   Kmd2 = Kmd_gas;
    umd1 = umd_oil;   umd2 = umd_gas;
end
% [Kf_wood] = wood(Kf(1),Kf(2),S_1);

Kmd1 =conj(Kmd1);   Kmd2 =conj(Kmd2);  
umd1 =conj(umd1);   umd2 =conj(umd2); 

%%  High and low frequency limit
[Ksat_water,usat_water] = Gassmann (Kmd_water,umd_water,Ks,por,Kf_water);
Q_sq_water = abs(imag((Ksat_water + 4./3.*usat_water))./real((Ksat_water + 4./3.*usat_water)));

[Ksat_gas,usat_gas] = Gassmann (Kmd_gas,umd_gas,Ks,por,Kf_gas);
Q_sq_gas = abs(imag((Ksat_gas + 4./3.*usat_gas)).*10./real((Ksat_gas + 4./3.*usat_gas)));

npress = [5];          %  Pressure（MPa）
np = npress.*5;      %  The serial number corresponding to the pressure
S1 = 0.8;       %  Given saturation
perm = 1e-15;   %  Given permeability
S2 = 1 - S1;
rho = rho_s*(1-por) + por*S1*rho_f(1) + por*S2*rho_f(2);
for j = 1:n2
    b = 0.05;     % Parameters related to the autocorrelation function
    [Hh1(:,j),Hw1(:,j),H_sc1(:,j),Vp_crm(:,j),Qp_crm(:,j)] = ...
        CRM(Kmd1(:,j)',Kmd2(:,j)',umd1(:,j)',umd2(:,j)',Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
    K_crm(:,j) = H_sc1(:,j) - 4.*umd(:,j)./3;

    [Hh2(:,j),Hw2(:,j),H_sc2(:,j),Vp_crm2(:,j),Qp_crm2(:,j)] = ...
        CRM(K_d(j),K_d(j),um_d(j),um_d(j),Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
end
Vp_GW1 = sqrt(Hw1./rho);
Vp_GW2 = sqrt(Hw2./rho);
Vp_GH1 = sqrt(Hh1./rho);
Vp_GH2 = sqrt(Hh2./rho);


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);                    
semilogx(FrequencyRange,Vp_crm(:,np(1))./1000,'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Vp_crm2(:,np(1))./1000,'-k','linewidth',2.0);
% semilogx(FrequencyRange,Vp_GH1(:,np(1))./1000,'--m','linewidth',2.0);
% semilogx(FrequencyRange,Vp_GW1(:,np(1))./1000,'--b','linewidth',2.0);
semilogx(FrequencyRange,Vp_GH2(np(1))./1000.*ones(1,length(FrequencyRange)),'--g','linewidth',2.0);
semilogx(FrequencyRange,Vp_GW2(np(1))./1000.*ones(1,length(FrequencyRange)),'--c','linewidth',2.0);
legend('Scaled','Patchy','Gassmann-Hill','Gassmann-Wood');
xlabel('Frequency[Hz]');
ylabel('{\itV}_P (km/s)');
xlim([1e-4 1e7]); 
ylim([3.600 4.400]);
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'ytick',[3.6 3.8 4 4.2 4.4],'fontsize',Fontsize1);
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
semilogx(FrequencyRange,Qp_crm(:,np(1)),'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Qp_crm2(:,np(1)),'-k','linewidth',2.0);
semilogx(FrequencyRange,Q_sq_water(:,np(1)),'-b','linewidth',2.0);
semilogx(FrequencyRange,Q_sq_gas(:,np(1)),'-m','linewidth',2.0);
legend('Scaled','Patchy','Squirt-Water','Squirt-Gas');
xlabel('Frequency[Hz]');
ylabel('1/Q_p');
xlim([1e-4 1e7]); 
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

%%  The influence of stress
% Given Pressure
npress = [1 10 20];          %  Pressure（MPa）
np = npress.*5;      %  The serial number corresponding to the pressure
S1 = 0.8;       %  Given saturation
perm = 1e-15;   %  Given permeability
S2 = 1 - S1;
rho = rho_s*(1-por) + por*S1*rho_f(1) + por*S2*rho_f(2);
for j = 1:n2
    b = 0.05;     % Parameters related to the autocorrelation function
    [Hh1(:,j),Hw1(:,j),H_sc1(:,j),Vp_crm(:,j),Qp_crm(:,j)] = ...
        CRM(Kmd1(:,j)',Kmd2(:,j)',umd1(:,j)',umd2(:,j)',Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
    K_crm(:,j) = H_sc1(:,j) - 4.*umd(:,j)./3;

    [Hh2(:,j),Hw2(:,j),H_sc2(:,j),Vp_crm2(:,j),Qp_crm2(:,j)] = ...
        CRM(K_d(j),K_d(j),um_d(j),um_d(j),Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);

    [Vp_js(:,j),Vs_js(:,j),Qp_js(:,j),Qs_js(:,j),K(:,j)]  =  Johnson_ball_new(Kmd1(:,j)',umd1(:,j)',Ks,rho_s,por,perm,Kf(1),Kf(2),vis_f(1),vis_f(2),S1,S2,rho_f(1),rho_f(2),w,0.1);
end
Vp_GW1 = sqrt(Hw1./rho);
Vp_GW2 = sqrt(Hw2./rho);
Vp_GH1 = sqrt(Hh1./rho);
Vp_GH2 = sqrt(Hh2./rho);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);                     
semilogx(FrequencyRange,Vp_crm(:,np(1))./1e3,'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Vp_crm(:,np(2))./1e3,'-b','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm(:,np(3))./1e3,'-k','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm2(:,np(1))./1e3,'--r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Vp_crm2(:,np(2))./1e3,'--b','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm2(:,np(3))./1e3,'--k','linewidth',2.0);
legend('press = 1MPa''press = 10MPa','press = 20MPa');
xlabel('Frequency[Hz]');
ylabel('{\itV}_P (km/s)');
xlim([1e-4 1e7]); 
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
semilogx(FrequencyRange,Qp_crm(:,np(1)),'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Qp_crm(:,np(2)),'-b','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm(:,np(3)),'-k','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm2(:,np(1)),'--r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Qp_crm2(:,np(2)),'--b','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm2(:,np(3)),'--k','linewidth',2.0);
legend('press = 1MPa','press = 10MPa','press = 20MPa');
xlabel('Frequency[Hz]');
% ylabel('Attenuation[1/Q_p]');
ylabel('1/Q_p');
xlim([1e-4 1e7]); 
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

%%   The influence of saturation
% Given Pressure
npress = 5;          %  Pressure（MPa）
np = npress.*5;      %  The serial number corresponding to the pressure（MPa）
perm = 1e-15;        %  Given permeability
for i = 1:length(S_1)
    S1 = S_1(i);
    S2 = 1 - S1;
    rho = rho_s*(1-por) + por*S1*rho_f(1) + por*S2*rho_f(2);
    for j = 1:n2
        b = 0.05;     % Parameters related to the autocorrelation function
        [Hh1(:,j),Hw1(:,j),H_sc1(:,j),Vp_crm(:,j),Qp_crm(:,j)] = ...
            CRM(Kmd1(:,j)',Kmd2(:,j)',umd1(:,j)',umd2(:,j)',Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
        K_crm(:,j) = H_sc1(:,j) - 4.*umd(:,j)./3;
    end
    Vp_crm_Sw(i,:) = Vp_crm(:,np);
    Qp_crm_Sw(i,:) = Qp_crm(:,np);
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);                    
semilogx(FrequencyRange,Vp_crm_Sw(1,:)./1000,'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Vp_crm_Sw(2,:)./1000,'-b','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_Sw(3,:)./1000,'-k','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_Sw(4,:)./1000,'-m','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_Sw(5,:)./1000,'-g','linewidth',2.0);
legend('Sw = 0.1','Sw = 0.3','Sw = 0.5','Sw = 0.7','Sw = 0.9');
xlabel('Frequency[Hz]');
ylabel('{\itV}_P (km/s)');
xlim([1e-4 1e7]); 
ylim([3.6 4.4]);
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
semilogx(FrequencyRange,Qp_crm_Sw(1,:),'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Qp_crm_Sw(2,:),'-b','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_Sw(3,:),'-k','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_Sw(4,:),'-m','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_Sw(5,:),'-g','linewidth',2.0);
legend('Sw = 0.1','Sw = 0.3','Sw = 0.5','Sw = 0.7','Sw = 0.9');
xlabel('Frequency[Hz]');
ylabel('1/Q_p');
xlim([1e-4 1e7]); 
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

% save Vp_Sw Vp_crm_Sw;
% save Qp_Sw Qp_crm_Sw;

%%  The influence of permeability
% Given Pressure
npress = 5;          %  Pressure（MPa）
np = npress.*5;      %  The serial number corresponding to the pressure（MPa）
S1 = 0.8;            %  Given saturation
S2 = 1 - S1;
rho = rho_s*(1-por) + por*S1*rho_f(1) + por*S2*rho_f(2);
for i = 1:length(perm0)
    perm = perm0(i);
    for j = 1:n2
        b = 0.05;     % Parameters related to the autocorrelation function
        [Hh1(:,j),Hw1(:,j),H_sc1(:,j),Vp_crm(:,j),Qp_crm(:,j)] = ...
            CRM(Kmd1(:,j)',Kmd2(:,j)',umd1(:,j)',umd2(:,j)',Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
        K_crm(:,j) = H_sc1(:,j) - 4.*umd(:,j)./3;
    end
    Vp_crm_perm(i,:) = Vp_crm(:,np);
    Qp_crm_perm(i,:) = Qp_crm(:,np);
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);                   
semilogx(FrequencyRange,Vp_crm_perm(1,:)./1000,'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Vp_crm_perm(2,:)./1000,'-b','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_perm(3,:)./1000,'-k','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_perm(4,:)./1000,'-m','linewidth',2.0);
semilogx(FrequencyRange,Vp_crm_perm(5,:)./1000,'-g','linewidth',2.0);
legend('perm = 0.1mD','perm = 1mD','perm = 10mD','perm = 100mD','perm = 1D');
xlabel('Frequency[Hz]');
ylabel('{\itV}_P (km/s)');
xlim([1e-4 1e7]); 
ylim([3.6 4.4]);
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
semilogx(FrequencyRange,Qp_crm_perm(1,:),'-r','linewidth',2.0);hold on;
semilogx(FrequencyRange,Qp_crm_perm(2,:),'-b','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_perm(3,:),'-k','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_perm(4,:),'-m','linewidth',2.0);
semilogx(FrequencyRange,Qp_crm_perm(5,:),'-g','linewidth',2.0);
legend('perm = 0.1mD','perm = 1mD','perm = 10mD','perm = 100mD','perm = 1D');
xlabel('Frequency[Hz]');
ylabel('1/Q_p');
xlim([1e-4 1e7]); 
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca,'xtick',[10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4) 10^(6)],'fontsize',Fontsize1);
set(gca, 'FontSize',Fontsize1);
grid on; grid minor;

% save Vp_perm Vp_crm_perm;
% save Qp_perm Qp_crm_perm;











