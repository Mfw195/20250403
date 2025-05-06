
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----             双尺度模型与实验数据标定        (部分饱和甘油)
clc;  clear all;
%%%%%%%%%%%%%%%%%%%%%%%     将 改进的骨架模量带入斑块模型中 并分析其敏感性  %%%%%%%%%%%%%%%%%%%%
Fontsize1 = 14;
load ./model_data./Kmd_gas.mat;
load ./model_data./umd_gas.mat;
Kmd_gas = Kmd;
umd_gas = umd;
load ./model_data./kmd_water.mat;
load ./model_data./umd_water.mat;
Kmd_water = Kmd;
umd_water = umd;
Q_md_water = abs(imag((Kmd_water + 4./3.*umd_water))./real((Kmd_water + 4./3.*umd_water)));
load ./model_data./kmd_oil.mat;
load ./model_data./umd_oil.mat;
Kmd_oil = Kmd;
umd_oil = umd;

load ./model_data./K_dry.mat;
load ./model_data./u_dry.mat;

FrequencyRange1=[10.^(-4:0.01:-0.01) 1:1:200 10.^(2.4:0.01:8)];
w=2*pi*FrequencyRange1;

Ks = 39e9;
us = 35e9;
rho_s = 2650;
por = 0.12;   
[n1,n2] = size(Kmd_oil);
%%%%%%%%%%%%%%% Parameters of fluids ----------------------------------------------------
Kf_oil     = 4.66*10^(9);    % 甘油
rho_oil = 1261;            % (Kg/m^3)
vis_oil    = 900*10^(-3);      % (pa.s) 1pa.s=1000mpa.s=1000cP

Kf_gas  = 0.01*10^(9);              % gas （衰减幅值反比，峰值频率正比）
rho_gas = 50;
vis_gas = 1*10^(-5);          

Kf_water     = 2.25*10^(9);        % Water
rho_water = 1000;                
vis_water    = 1*10^(-3);
FrequencyRange = FrequencyRange1.*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 3;                         %  选择1--水+气，2--水+油， 3--油+气
S_1 = [0.1 0.3 0.5 0.7 0.9];   %  S_water : 流体1的饱和度
perm0 = [1e-17 1e-16 1e-15 1e-14 1e-13] ;

if n ==1
    Kf = [Kf_water Kf_gas];                 % 各流体体积模量
    rho_f = [rho_water rho_gas];            % 各流体的密度
    vis_f = [vis_water vis_gas];            % 各流体的粘度
    Kmd1 = Kmd_water;   Kmd2 = Kmd_gas;
    umd1 = umd_water;   umd2 = umd_gas;
elseif n == 2
    Kf = [Kf_water Kf_oil];                 % 各流体体积模量
    rho_f = [rho_water rho_oil];            % 各流体的密度
    vis_f = [vis_water vis_oil];            % 各流体的粘度
    Kmd1 = Kmd_water;   Kmd2 = Kmd_oil;
    umd1 = umd_water;   umd2 = umd_oil;
else
    Kf = [Kf_oil Kf_gas];                 % 各流体体积模量
    rho_f = [rho_oil rho_gas];            % 各流体的密度
    vis_f = [vis_oil vis_gas];            % 各流体的粘度
    Kmd1 = Kmd_oil;   Kmd2 = Kmd_gas;
    umd1 = umd_oil;   umd2 = umd_gas;
end
% [Kf_wood] = wood(Kf(1),Kf(2),S_1);

Kmd1 =conj(Kmd1);   Kmd2 =conj(Kmd2);  
umd1 =conj(umd1);   umd2 =conj(umd2); 

% 选择压力
npress = [0.2 3 5 8.8 12 15 20];          %  压力（MPa）
np = npress.*5;      %  压力对应的序号（MPa）

S1 = 0.95;       %  给定饱和度
perm = 1e-15;   %  给定渗透率
S2 = 1 - S1;
rho = rho_s*(1-por) + por*S1*rho_f(1) + por*S2*rho_f(2);
for j = 1:n2
    b = 0.004;     % 与自相关函数有关的参数   与特征频率反比
    [Hh1(:,j),Hw1(:,j),H_sc1(:,j),Vp_crm1(:,j),Qp_crm(:,j)] = ...
        CRM(Kmd1(:,j)',Kmd2(:,j)',umd1(:,j)',umd2(:,j)',Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
    K_crm(:,j) = H_sc1(:,j) - 4.*umd(:,j)./3;

    Vp_crm(:,j) = fillmissing(Vp_crm1(:,j),'previous');

    [Hh2(:,j),Hw2(:,j),H_sc2(:,j),Vp_crm2(:,j),Qp_crm2(:,j)] = ...
        CRM(K_d(j),K_d(j),um_d(j),um_d(j),Ks,por,S1,S2,Kf(1),Kf(2),vis_f(1),vis_f(2),rho,w,perm,b,1);
end



%%%%%%%%%%%%%%%%%%%%%%%%  读取实验数据 --- 甘油
%%%%%%%%%%%%%%%%%%   =====================================================
CC = colormap(jet(8)); Fontsize1=14;

load ./Experimental_data/Low_frequency_data/Vp_low_glycerol.mat;
load ./Experimental_data/Low_frequency_data/Qp_low_glycerol.mat;
load ./Experimental_data/Low_frequency_data/F.mat;

figure(1)
for i = 1:7
    semilogx(FrequencyRange,Vp_crm(:,np(i))./1e3,'Color',CC(i,:),'linewidth',2.0);  hold on;
end
for i = 1:7
    semilogx(freq,Vp_W(:,i)./1000,'o','LineWidth',1.3,'Color',CC(i,:)); hold on;
end
%ylim([2400 4500]);
xlabel('Frequency (Hz)','Fontsize',Fontsize1);
ylabel('{\itV}_P (km/s)','Fontsize',Fontsize1);
set(gca,'xtick',[10^(-2),10^(-1),1,10,10^2,10^3,10^4,10^5,10^6],'Fontsize',Fontsize1);
%set(gca,'ytick',[20 25 30 35 40 45 50 55 60 65 70 75 80],'Fontsize',Fontsize1);
xlim([10^(-1) 10^7]);
ylim([3 5]);
title('(b)','FontName','Times New Roman','Fontsize',Fontsize1);
legend('0 MPa','3MPa','5 MPa','8MPa','10 MPa','15 MPa','20 MPa','FontName','Times New Roman','Fontsize',Fontsize1);
box on;     grid on;  grid minor;

figure(2)
for i = 1:2:7
    semilogx(FrequencyRange,Qp_crm(:,np(i)),'Color',CC(i,:),'linewidth',2.0);  hold on;
end
for i = 1:2:7
    semilogx(freq(1:34),QP_W(:,i),'o','LineWidth',1.3,'Color',CC(i,:)); hold on;
end
xlabel('Frequency (Hz)','Fontsize',Fontsize1);
ylabel('1/Q_p','Fontsize',Fontsize1);
set(gca,'xtick',[1,10,10^2,10^3,10^4,10^5,10^6],'Fontsize',Fontsize1);
%set(gca,'ytick',[20 25 30 35 40 45 50 55 60 65 70 75 80],'Fontsize',Fontsize1);
xlim([10^0 10^4]);
ylim([0 0.12]);
title('(d)','FontName','Times New Roman','Fontsize',Fontsize1);
legend('0MPa','5 MPa','10 MPa','20 MPa','FontName','Times New Roman','Fontsize',Fontsize1);
box on;     grid on;  grid minor;

