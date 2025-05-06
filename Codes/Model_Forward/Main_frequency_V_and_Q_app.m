clc;  clear all;
%%%%%%%%%%%%%%%%%%%%%%%     将 改进的骨架模量带入斑块模型中 并分析其敏感性  %%%%%%%%%%%%%%%%%%%%
tic;
red = [1 0 0];  green = [0 1 0];  blue = [0 0 1];   magenta = [1 0 1];   black = [0 0 0];   
CC = [red;blue;black;magenta;green];

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

FrequencyRange=[10.^(-4:0.01:-0.01) 1:1:200 10.^(2.4:0.01:8)];
w=2*pi*FrequencyRange;

Ks = 49e9;
us = 35e9;
rho_s = 2650;
por = 0.12;   
[n1,n2] = size(Kmd_oil);
%%%%%%%%%%%%%%% Parameters of fluids ----------------------------------------------------
Kf_oil     = 4.66*10^(9);    % 甘油
rho_oil = 1261;            % (Kg/m^3)
vis_oil    = 1400*10^(-3);      % (pa.s) 1pa.s=1000mpa.s=1000cP

Kf_gas  = 0.01*10^(9);              % gas （衰减幅值反比，峰值频率正比）
rho_gas = 50;
vis_gas = 1*10^(-5);          

Kf_water     = 2.25*10^(9);        % Water
rho_water = 1000;                
vis_water    = 1*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;                         %  选择1--水+气，2--水+油， 3--油+气
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

%%   饱和度和渗透率共同作用
% 选择压力
npress = 5;          %  压力（MPa）
np = npress.*5;      %  压力对应的序号（MPa）

S_1 = 0.01:0.01:0.99;   %  S_water : 流体1的饱和度
perm0 = 10.^(-16:0.01:-12) ;

S_2 = 1 - S_1;
rho = rho_s*(1-por) + por.*S_1.*rho_f(1) + por.*S_2.*rho_f(2);
b = 0.05;     % 与自相关函数有关的参数   与特征频率反比
ze = zeros(length(S_1),length(perm0),length(w));
Vp_crm_1 = ze;
Qp_crm_1 = ze;
for i = 1:length(S_1)
    for j = 1:length(perm0)
  
        [Hh1,Hw1,H_sc1,Vp_crm,Qp_crm] = ...
            CRM(Kmd1(:,np)',Kmd2(:,np)',umd1(:,np)',umd2(:,np)',Ks,por,S_1(i),S_2(i),Kf(1),Kf(2),vis_f(1),vis_f(2),rho(i),w,perm0(j),b,1);
        Vp_crm_1(i,j,:) = Vp_crm;
        Qp_crm_1(i,j,:) = Qp_crm;
    end
end

% save('./model_data./Vp_Sw_perm.mat','Vp_crm_1') ;
% save('./model_data./Qp_crm_perm.mat','Qp_crm_1') ;

values = find(ismember(FrequencyRange,30));
Vp_sq_crm_30Hz = abs(Vp_crm_1(:,:,values));
Qp_sq_crm_30Hz = Qp_crm_1(:,:,values);
CC1=colormap(jet(5));    CC2=colormap(hot(4));   CC3=colormap(cool(5));   CC4=colormap(spring(7));

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
imagesc(perm0.*1e15,S_1.*1e2,Vp_sq_crm_30Hz);   hold on;
% yline(nb,'r',LineWidth=1.5);
ylabel('{\itS}_W(%)','FontName','Times New Roman');
xlabel('Permeability(mD)','FontName','Times New Roman');
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
axis xy;
set(gca,'ytick',[10 30 50 70 90],'FontSize',14);
set(gca,'xScale','log','xtick',[1e-1 1 10 100 10^3],'FontSize',14);     %  Vp分布图
grid on
ylim([10 90]);
xlim([0.1 1000]);
colormap('hot');
colorbar;
my_handle=colorbar;
my_handle.Title.String = '{\itV}_P (30Hz)[m/s]';
my_handle.Title.FontSize = 12;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]); 
imagesc(perm0.*1e15,S_1.*1e2,Qp_sq_crm_30Hz);      hold on;
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
ylabel('{\itS}_W(%)','FontName','Times New Roman');
xlabel('Permeability(mD)','FontName','Times New Roman');
% title('频率30Hz时的1/Q');
axis xy;
set(gca,'ytick',[10 30 50 70 90],'FontSize',14);
set(gca,'xScale','log','xtick',[1e-1 1 10 100 1000],'FontSize',14);     %  Q值分布图
grid on
ylim([10 90]);
xlim([0.1 1000]);
colormap('hot');
colorbar;
my_handle=colorbar;
my_handle.Title.String = '1/{\itQ}_P (30Hz)';
my_handle.Title.FontSize = 12;

Sw = S_1.*100;
n_Sw = [90];
N_Sw = find(ismember(Sw,n_Sw));   
figure
for i = 1:length(n_Sw)
semilogx(perm0.*1e15,Vp_sq_crm_30Hz(N_Sw(i),:),'color',CC(i,:),LineWidth=1.5);  hold on;
end
xlim([0.1 1000]);
set(gca,'XScale','log','xtick',[0.1 1 10 100 1000],'FontSize',14);
xlabel('Permeability (mD)');
ylabel('Velocity (m/s)');
legend('Sw = 80%');
grid on;   grid minor;  box on;
set(gca, 'FontSize',Fontsize1);
set(gcf, 'Position', [10 220 500 300]);

figure
for i = 1:length(n_Sw)
semilogx(perm0.*1e15,Qp_sq_crm_30Hz(N_Sw(i),:),'color',CC(i,:),LineWidth=1.5);  hold on;
end
xlabel('Permeability (mD)');
ylabel('1/Q');
xlim([0.1 1000]);
ylim([0.01 0.045]);
set(gca,'XScale','log','xtick',[0.1 1 10 100 1000],'FontSize',14);
legend('Sw = 80%');
grid on;   grid minor;  box on;
set(gca, 'FontSize',Fontsize1);
set(gcf, 'Position', [10 220 500 300]);

perm = perm0.*1e18;  %  uD
n_perm = [1e-15].*1e18;
N_perm = find(ismember(perm,n_perm));   

figure
for i = 1:length(n_perm)
plot(Sw,Vp_sq_crm_30Hz(:,N_perm(i)),'color',CC(i,:),LineWidth=1.5);  hold on;
end
xlabel('Sw (%)');
ylabel('P-wave velocity (m/s)');
legend('Permeability = 1mD');
grid on;   grid minor;  box on;
set(gca, 'FontSize',Fontsize1);
set(gcf, 'Position', [10 220 500 300]);

figure
for i = 1:length(n_perm)
plot(Sw,Qp_sq_crm_30Hz(:,N_perm(i)),'color',CC(i,:),LineWidth=1.5);  hold on;
end
xlabel('Sw (%)');
ylabel('1/Q_p');
legend('Permeability = 1mD');
grid on;   grid minor;  box on;
set(gca, 'FontSize',Fontsize1);
set(gcf, 'Position', [10 220 500 300]);









