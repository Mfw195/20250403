clc;    clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================     根据超声数据反演得到纵横比分布    ========================================================================
Fontsize = 14;              % For display
AspectRatio = 0.618; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_frame = 2820;                                  % Density of solid grains: kg/m^3（高低速度频极限，反比）
np=7;  CC=colormap(jet(np));
%%%%%%% 岩石骨架速度随压力的变化 ============================================
p = [100000	1000000	3000000	5000000	8000000	10000000	15000000	20000000	...
    25000000	30000000	35000000	40000000	45000000	50000000];           % Effective pressure, pa
press=p/10^6;                % Effective pressure, Mpa

vp_dry = [2945.94594594595	3103.20284697509	3303.03030303030	3433.07086614173	3603.30578512397	3710.63829787234	3945.70135746606	...
    4093.89671361502	4172.24880382775	4295.56650246305	4360	4404.04040404040	4426.39593908630	4518.13471502591];
vs_dry = [2088.62275449102	2101.20481927711	2153.08641975309	2218.82951653944	2294.73684210526	2356.75675675676	2463.27683615819	...
    2557.18475073314	2618.61861861862	2674.84662576687	2716.51090342679	2742.13836477987	2768.25396825397	2803.85852090032];

err_Vp = [0.0936784713107621	0.130775529616314	0.0842219254301046	0.127928436884867	0.136513511862060	0.0676149074817211	0.0712515924158343	0.0851558712301278...
    0.107640372837455	0.114089908539822	0.0714163673647971	0.115193879313703	0.0799634626379405	0.0906038380129680];
err_Vs = [0.0720397601291064	0.0890395910787588	0.0467792876313422	0.0817697618387378	0.0750192627604094	0.0677542252140319	0.0631439980506347	0.0583942037133012...
    0.0423666648295219	0.0855674792993412	0.0665166698837150	0.0550244460590489	0.0402269110615336	0.0970339089307946];

Fontsize1 = 14;
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(press, vp_dry./1000,'--bd','MarkerFaceColor','b','MarkerSize',7,LineWidth=1.2);  hold on;
plot(press, vs_dry./1000,'--r^','MarkerFaceColor','r','MarkerSize',7,LineWidth=1.2);  hold on;

errorbar(press,vp_dry./1000,err_Vp,'--bd'); hold on;
errorbar(press,vs_dry./1000,err_Vs,'--r^');

xlabel('Effective pressure {\itP_{eff}} (MPa)','Fontsize',Fontsize1);
ylabel('Velocity (km/s)','Fontsize',Fontsize1);
set(gca,'xtick',[0 10 20 30 40 50],'Fontsize',Fontsize1);
set(gca,'ytick',[2 2.5 3 3.5 4 4.5 5],'Fontsize',Fontsize1);
legend('P-wave','S-wave');
ylim([2 5]);
xlim([-2 52]);
title('(a)','FontSize',Fontsize1,'FontName','Times New Roman');
grid on;   box on;  grid minor;

%%%%%%%% Moduli calculations based on parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%------------------超声实验测得的岩石弹性模量---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_dry_c=rho_frame;                       % Density of solid grains
um_dry_c = rho_dry_c*vs_dry.^2;               % Shear modulus
lamda_dry_c = rho_dry_c*(vp_dry.^2-2*vs_dry.^2);     % Lamda
Kdry_c  = lamda_dry_c+2*um_dry_c/3;                % Bulk modulus
po_dry_c = 0.5*lamda_dry_c./(lamda_dry_c+um_dry_c);    % The Poisson's ratio

%%% espilon, 岩石中所有未闭合软孔隙的累积裂缝密度
espilon1=(Kdry_c(end)./Kdry_c-1)*9*(1-2*po_dry_c(end))/16/(1-po_dry_c(end)^2);                % Crack density from the bulk modulus, in Eq.(5) of 2015_CG
espilon2=(um_dry_c(end)./um_dry_c-1)*45*(2-po_dry_c(end))/32/(1-po_dry_c(end))/(5-po_dry_c(end));  % Crack density from the shear modulus, in Eq.(5) of 2015_CG

z1=log(espilon1(1:end-1));       %%% Fitting the curves  （拟合）
z2=log(espilon2(1:end-1));       %
c1=polyfit(p(1:end-1),z1,1);
c2=polyfit(p(1:end-1),z2,1);     %（拟合裂缝密度与压力的关系）

dp=0.2*10^6;           % Pressure interval, 0.2*10^6;
pp=0:dp:p(end);        % Pressure, MPa
press_pp=pp/10^6;      % Effective pressure, MPa

ie1=exp(c1(1)*pp+c1(2));         %%% 不同压力下的裂缝密度
ie2=exp(c2(1)*pp+c2(2));         % Fitting crack density, Eq.(7) of 2015_CG

ie = (ie1+ie2)/2;                %%% Mean of bulk and shear moduli

% K_d & um_d: the dry effective modulus; %%% 整个压力范围内，计算拟合后的干骨架模量
K_d  = Kdry_c(end)./(1+16*(1-po_dry_c(end)^2)*ie1/9/(1-2*po_dry_c(end)));           % Kd in Eq.(5) in 2015_CG, Eq.(8) in 2018_JGE;
um_d = um_dry_c(end)./(1+32*(1-po_dry_c(end))*(5-po_dry_c(end))*ie2/45/(2-po_dry_c(end)));
%%% Kdry_c(end):   Kdry_ciff of Eq.(5) in 2015_CG
%%% po_ra(end):  The Poisson's ratio of the background matrix with only hard pores

lamda_d=K_d-2*um_d/3;
E_d=3*K_d.*um_d./(lamda_d+um_d);
po_d = 0.5*lamda_d./(lamda_d+um_d);        %拟合后干岩石模量

lamda_st = Kdry_c(end)-2*um_dry_c(end)/3;                   % High-freq limit；%%% 计算最大压力下的模量，软孔隙接近完全闭合
E_st=3*Kdry_c(end).*um_dry_c(end)./(lamda_st+um_dry_c(end));
po_st = 0.5*lamda_st./(lamda_st+um_dry_c(end));
K_st = Kdry_c(end);
um_st = um_dry_c(end);

alpha_i = 4*(1-po_st.^2).*pp./(E_st*pi);               % Eq.(10) in 2015_CG, Eq.(12) in 2018_JGE
%%% 各个纵横比随压力的变化（含裂隙分布岩石中未闭合软孔隙的最小初始纵横比

%%% Aspect ratio Variations with pressure: Eq.(13) in 2018_JGE
alpha_p=zeros(length(pp),length(pp));      %% 所有软孔隙的纵横比值与压力的关系：

iei=[0 ie(1:end-1)-ie(2:end)];       % 各个纵横比(or different pressure)对应的裂隙密度
iek=[0 ie1(1:end-1)-ie1(2:end)];     % 各个纵横比对应的裂隙密度
ieu=[0 ie2(1:end-1)-ie2(2:end)];     % 各个纵横比对应的裂隙密度

for ix=1:length(pp)
    if( ix==1 )
        alpha_p(:,ix) = alpha_i-4*(1-po_st^2)*0*dp/(3*pi*K_st*(1-2*po_st));     % Eq.(11) in 2015_CG, Eq.(13) in 2018_JGE
    elseif (ix>1)
        alpha_p(:,ix) = alpha_p(:,ix-1)-4*(1-po_st^2)*dp/(3*pi*K_st*(1-2*po_st));
    end
end

iei1=iei'*ones(1,length(pp));      %% 累积裂隙密度分布图; 算法二
iek1=iek'*ones(1,length(pp));
ieu1=ieu'*ones(1,length(pp));
phy_k=4/3*pi*alpha_p.*iei1;            % Line-15 in P-3394 in 2015_CG

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);               % 软孔隙累积裂隙密度
plot(press,espilon1,'dk','MarkerFaceColor','k','MarkerSize',7);    hold on;
plot(press,espilon2,'^b','MarkerFaceColor','b','MarkerSize',7);
plot(press_pp,ie,'-r','linewidth',2.0);
legend('From bulk modulus','From shear modulus','Fitting');
set(gca,'xtick',[0 10 20 30 40 50 60],'fontsize',Fontsize);
set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1.0],'fontsize',Fontsize);
xlabel('Effective pressure {\itP_{eff}} (MPa)','fontsize',Fontsize);
ylabel('Crack density','fontsize',Fontsize);
xlim([-2 52]);   ylim([-0.05 1]);
title('(b)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca, 'FontSize',Fontsize1) ;
grid on;  grid minor;  box on;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(alpha_p(1:5:end,1),  phy_k(1:5:end,1).*1e2,'-k','linewidth',2.0);   hold on
plot(alpha_p(1:5:end,5e6/dp+1),phy_k(1:5:end,5e6/dp+1).*1e2,'-m','linewidth',2.0);
plot(alpha_p(:,10e6/dp+1),phy_k(:,10e6/dp+1).*1e2,'-b','linewidth',2.0);
plot(alpha_p(:,20e6/dp+1),phy_k(:,20e6/dp+1).*1e2,'-r','linewidth',2.0);
plot(alpha_p(:,40e6/dp+1),phy_k(:,40e6/dp+1).*1e2,'-g','linewidth',2.0);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2].*1e-3,'fontsize',Fontsize);
set(gca,'ytick',[0 1.0 2.0 3 4 5].*1e-4,'fontsize',Fontsize);
xlabel('Aspect ratio (\alpha)','fontsize',Fontsize);      % 软孔隙纵横比
ylabel('Crack porosity (%) ','fontsize',Fontsize);       % 软孔隙度
xlim([0 1.3e-3]);       ylim([0 5e-4]);
title('(c)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca, 'FontSize',Fontsize1) ;
legend('0 MPa','5 MPa','10 MPa','20 MPa','40 MPa');
grid on;  grid minor;  box on;

nn1=1;nn2=1;                        % 寻找每一压力下软孔隙度和纵横比非零值的起点
mn=zeros(1,length(pp));
for np=1:length(pp)
    for ix=1:length(pp)
        if phy_k(ix,np)<=1e-8
            nn1=ix+1;
        end
        if alpha_p(ix,np)<=1e-8
            nn2=ix+1;
        end
    end
    mn(np)=max([nn1 nn2]);
end

for np=1:length(pp)
    for ix=mn(np):length(pp)
        iei1(ix,np)=iei1(ix-1,np)+iei1(ix,np);
        iek1(ix,np)=iek1(ix-1,np)+iek1(ix,np);
        ieu1(ix,np)=ieu1(ix-1,np)+ieu1(ix,np);     %累计裂隙密度
    end
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(alpha_p(1:5:end,1),  iei1(1:5:end,2),'-k','linewidth',2.0);hold on
plot(alpha_p(:,5e6/dp+1),iei1(:,5e6/dp+1),'-m','linewidth',2.0);
plot(alpha_p(:,10e6/dp+1),iei1(:,10e6/dp+1),'-b','linewidth',2.0);
plot(alpha_p(:,20e6/dp+1),iei1(:,20e6/dp+1),'-r','linewidth',2.0);
plot(alpha_p(:,40e6/dp+1),iei1(:,40e6/dp+1),'-g','linewidth',2.0);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2]./1e3,'fontsize',Fontsize);
set(gca,'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7],'fontsize',Fontsize);
xlabel('Aspect ratio (\alpha/10^3)','fontsize',Fontsize);      % 软孔隙纵横比
ylabel('Cumulative crack density');             % 未闭合软孔隙累计裂隙密度
xlim([0 1.3e-3]);        ylim([0 0.7]);
title('(d)','FontSize',Fontsize1,'FontName','Times New Roman');
set(gca, 'FontSize',Fontsize1) ;
legend('0 MPa','5MPa','10 MPa','20 MPa','40 MPa');
grid on;  grid minor;  box on;

% save alpha_distribution alpha_p;
% save phy_c_distribution phy_k;
% save press press_pp;
% save K_dry K_d;
% save u_dry um_d;