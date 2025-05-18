clc;    clear all;   close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----The influence of saturation on the forward simulation results —— AVA effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%============================================
Fontsize = 14;              % For display
AspectRatio = 0.618;
CC1=colormap(jet(7));    CC2=colormap(hot(7));   CC3=colormap(cool(7));   CC4=colormap(spring(7));
red = [1 0 0];  green = [0 1 0];  blue = [0 0 1];   magenta = [1 0 1];   black = [0 0 0];   
CC = [red;blue;black;magenta;green];

load ./model_data./Vp_Sw.mat;
load ./model_data./Qp_Sw.mat;
load ./model_data./Vp_perm.mat;
load ./model_data./Qp_perm.mat;

FrequencyRange=[10.^(-4:0.01:-0.01) 1:1:200 10.^(2.4:0.01:8)];
w=2*pi*FrequencyRange;

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
%%
% ===============================   Forward simulation      ===================================================  
V_forward = conj(Vp_crm_Sw);
[m,n] = size(V_forward);

Frequency1 = 1:1:100;
values = find(ismember(FrequencyRange,Frequency1));
for i =1:length(Frequency1)
    Vp_forward1(:,i) = V_forward(:,values(i));          
end
for i = 1:m
    Vp_0Hz(i) = min(V_forward(i,:));        
    Vp_forward(i,:) = [Vp_0Hz(i),Vp_forward1(i,:)];
end

Frequency = [0,Frequency1];

%  The velocities and densities of each elastic layer
% Vp1 = 4;
% Vp2 = 4.5;
% Vp4 = 4.1;
% Vs1 = 2.3;
% Vs2 = 2.5;
% Vs4 = 2.5;
% rho1 = 2.3;
% rho2 = 2.7;
% rho4 = 2.7;
Vp1 = 3.9;
Vp2 = 4.1;
Vp4 = 4.1;
Vs1 = 2.1;
Vs2 = 2.3;
Vs4 = 2.3;
rho1 = 2.2;
rho2 = 2.5;
rho4 = 2.5;
for i = 1:length(Frequency)
    Vp_1(i,:) = [Vp1 ,Vp2 , Vp_forward(1,i)./1000 ,Vp4];             %  The P-wave velocity of each layer(km/s)
    Vp_2(i,:) = [Vp1 ,Vp2 , Vp_forward(2,i)./1000 ,Vp4];             
    Vp_3(i,:) = [Vp1 ,Vp2 , Vp_forward(3,i)./1000 ,Vp4];             
    Vp_4(i,:) = [Vp1 ,Vp2 , Vp_forward(4,i)./1000 ,Vp4];            
    Vp_5(i,:) = [Vp1 ,Vp2 , Vp_forward(5,i)./1000 ,Vp4];             
end

Vs_rock = 2200.*ones(1,5);
Vs_1 = [Vs1 Vs2 Vs_rock(1)./1000 Vs4];             %  The S-wave velocity of each layer(km/s)
Vs_2 = [Vs1 Vs2 Vs_rock(2)./1000 Vs4];             
Vs_3 = [Vs1 Vs2 Vs_rock(3)./1000 Vs4];             
Vs_4 = [Vs1 Vs2 Vs_rock(4)./1000 Vs4];             
Vs_5 = [Vs1 Vs2 Vs_rock(5)./1000 Vs4];             

rho_s = 2650;
por = 0.12;   
S1 = [0.1 0.3 0.5 0.7 0.9];   
rho_f = [rho_water rho_gas];            
rho = rho_s.*(1-por) + por.*S1.*rho_f(1) + por.*(1 - S1).*rho_f(2);
rho_1 = [rho1 rho2 rho(1)./1000 rho4];          %  The density of each layer(g/cm^3)
rho_2 = [rho1 rho2 rho(2)./1000 rho4];          
rho_3 = [rho1 rho2 rho(3)./1000 rho4];          
rho_4 = [rho1 rho2 rho(4)./1000 rho4];          
rho_5 = [rho1 rho2 rho(5)./1000 rho4];         

thick = [0.2 0.2];          %  The thickness of each layer(km)
nL = length(thick);               %  Number of layers

d_ang=1;                                      % Angle interval, degree
ang_total = 30;                               % Maximum angle, degree
ang_inc=(0:d_ang:ang_total);                  % Incident angle: 1
theta =(0:d_ang:ang_total)*pi/180;            % Incident angle: 2

% -----------  The first saturation  ----------------------
for i = 1:length(theta)
    for j = 1:length(Frequency)
       [Rm] = ref_multi (Vp_1(1),Vs_1(1),rho_1(1),Vp_1(j,2:end),Vs_1(2:end),rho_1(2:end),thick,theta(i),Frequency(j));
       Rpp_1(i,j) = Rm(1,1);
    end
end
% -----------  The second saturation  ----------------------
for i = 1:length(theta)
    for j = 1:length(Frequency)
       [Rm] = ref_multi (Vp_2(1),Vs_2(1),rho_2(1),Vp_2(j,2:end),Vs_2(2:end),rho_2(2:end),thick,theta(i),Frequency(j));
       Rpp_2(i,j) = Rm(1,1);
    end
end
% -----------  The third saturation  ----------------------
for i = 1:length(theta)
    for j = 1:length(Frequency)
       [Rm] = ref_multi (Vp_3(1),Vs_3(1),rho_3(1),Vp_3(j,2:end),Vs_3(2:end),rho_3(2:end),thick,theta(i),Frequency(j));
       Rpp_3(i,j) = Rm(1,1);
    end
end
% -----------  The fourth saturation  ----------------------
for i = 1:length(theta)
    for j = 1:length(Frequency)
       [Rm] = ref_multi (Vp_4(1),Vs_4(1),rho_4(1),Vp_4(j,2:end),Vs_4(2:end),rho_4(2:end),thick,theta(i),Frequency(j));
       Rpp_4(i,j) = Rm(1,1);
    end
end
% -----------  The fifth saturation  ----------------------
for i = 1:length(theta)
    for j = 1:length(Frequency)
       [Rm] = ref_multi (Vp_5(1),Vs_5(1),rho_5(1),Vp_5(j,2:end),Vs_5(2:end),rho_5(2:end),thick,theta(i),Frequency(j));
       Rpp_5(i,j) = Rm(1,1);
    end
end

% ----     Frequency domain wavelet   -------------------
f_p = 30;           %  dominant frequency
tm = 0.05;          %  wavelet delay
N = 10000;    
dt = 1./N;
[WW_wavelet]= ricker2(f_p,Frequency,tm);
plot(Frequency,real(WW_wavelet));
xlim([0 100]);

%%%%% ----------------- Seismic synthetics ----============================
for i=1:length(theta)
    for j=1:length(Frequency) 
        syna_pp_r_1(j,i) = (Rpp_1(i,j))*WW_wavelet(j)  ; 
        syna_pp_r_2(j,i) = (Rpp_2(i,j))*WW_wavelet(j)  ; 
        syna_pp_r_3(j,i) = (Rpp_3(i,j))*WW_wavelet(j)  ;
        syna_pp_r_4(j,i) = (Rpp_4(i,j))*WW_wavelet(j)  ; 
        syna_pp_r_5(j,i) = (Rpp_5(i,j))*WW_wavelet(j)  ;
    end
end

for i = 1:length(theta)
    [t1,syna_p_r_1(:,i)] = FFT1(syna_pp_r_1(:,i),dt,-1);
    [t1,syna_p_r_2(:,i)] = FFT1(syna_pp_r_2(:,i),dt,-1);
    [t1,syna_p_r_3(:,i)] = FFT1(syna_pp_r_3(:,i),dt,-1);
    [t1,syna_p_r_4(:,i)] = FFT1(syna_pp_r_4(:,i),dt,-1);
    [t1,syna_p_r_5(:,i)] = FFT1(syna_pp_r_5(:,i),dt,-1);
end

tt1=(0:N-1)*dt*10^3;                 % Time for synthetics, ms

nd = [1 11 21 31];     %  Extract several channels to detect the amplitude differences
scale = 1;

Fontsize1 = 12;
%  ============  The first saturation =====================
% figure,
% for i = 1:length(nd)
%     plot(tt1,syna_p_r_1(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
% end
% ylim([-0.06 0.06]);
% set(gca,'ytick',[-0.04 0 0.04],'fontsize',Fontsize1);
% xlim([0 350])
% xlabel('Time (ms)');
% ylabel('amplitude');
% legend('0°','10°','20°','30°');
% grid on
% set(gcf, 'Position', [10 220 500 200]);

% %AVA
[amp_layer1,index1] = max(abs(syna_p_r_1(1:1000,:)),[],1);
[amp_layer2,index2] = max(abs(syna_p_r_1(1000:2000,:)),[],1);
[amp_layer3,index3] = max(abs(syna_p_r_1(2000:3000,:)),[],1);

% for i = 1:length(amp_layer1)
% syna_p_r_1(:,i) = syna_p_r_1(:,i)./max(abs(syna_p_r_1(:,i)));
% end
% for i = 1:length(index1)
%     amp_layer_1(i) = syna_p_r_1(index1(i),i);
%     amp_layer_2(i) = syna_p_r_1(1000+index2(i),i);
%     amp_layer_3(i) = syna_p_r_1(2000+index3(i),i);
% end
% amp_1ayar = [amp_layer_1 ; amp_layer_2 ; amp_layer_3 ];
% figure;
% for i = 1:3
%     scatter(ang_inc,amp_1ayar(i,:),'color',CC(i,:),LineWidth=1);hold on;
% end
% xlabel('Angle of incidence (°)');
% ylabel('Amplitude');
% % ylim([0 0.06]);
% % set(gca,'ytick',[0 0.02 0.04 0.06],'fontsize',Fontsize1);
% grid on
% set(gcf, 'Position', [10 220 500 200]);
time1 = tm.*1e3;
time2 = (thick(1)./Vp2.*1e3).*2 + time1;
time3 = (real(thick(2)./Vp_forward(1,30)).*1e6 + thick(1)./Vp2.*1e3).*2 + time1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.35 scrsz(4)*0.5]); hold on;
wigb(syna_p_r_1(:,1:2:end),scale,ang_inc(1:2:end),tt1);    % 反射率法
% wiggle(syna_p_r_1(:,1:2:end),tt1,ang_inc(1:2:end),1,'d');
yline(time1,'r--',LineWidth=2);
yline(time2,'r--',LineWidth=2);
yline(time3,'r--',LineWidth=2);
ylim([0 350]);
title('(a)','FontName','Times New Roman','Fontsize',Fontsize1);
xlabel('Incident angle (Degree)');
ylabel('Time (ms)');
% set(gcf, 'Position', [50 50 300 500]);
box on;

% ==============  The second saturation =====================
% figure,
% for i = 1:length(nd)
%     plot(tt1,syna_p_r_2(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
% end
% ylim([-0.2 0.3]);
% xlim([0 450])
% xlabel('Time (ms)');
% ylabel('amplitude');
% legend('0°','10°','20°','30°');
% grid on
% set(gcf, 'Position', [10 220 500 200]);
% 
% figure;hold on;
% wigb(syna_p_r_2,scale,ang_inc,tt1);    % 反射率法
% ylim([0 350]);
% xlabel('Incident angle (Degree)');
% ylabel('Time');

% =-============ The third saturation =====================
% figure,
% for i = 1:length(nd)
%     plot(tt1,syna_p_r_3(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
% end
% ylim([-0.06 0.06]);
% set(gca,'ytick',[-0.04 0 0.04],'fontsize',Fontsize1);
% xlim([0 350])
% xlabel('Time (ms)');
% ylabel('amplitude');
% legend('0°','10°','20°','30°');
% grid on
% set(gcf, 'Position', [10 220 500 200]);

[amp_layer1,index1] = max(abs(syna_p_r_3(1:1000,:)),[],1);
[amp_layer2,index2] = max(abs(syna_p_r_3(1000:2000,:)),[],1);
[amp_layer3,index3] = max(abs(syna_p_r_3(2000:3000,:)),[],1);

% for i = 1:length(amp_layer1)
% syna_p_r_3(:,i) = syna_p_r_3(:,i)./max(abs(syna_p_r_3(:,i)));
% end
% for i = 1:length(index1)
%     amp_layer_1(i) = syna_p_r_3(index1(i),i);
%     amp_layer_2(i) = syna_p_r_3(1000+index2(i),i);
%     amp_layer_3(i) = syna_p_r_3(2000+index3(i),i);
% end
% amp_1ayar = [amp_layer_1 ; amp_layer_2 ; amp_layer_3 ];
% figure;
% for i = 1:3
%     scatter(ang_inc,amp_1ayar(i,:),'color',CC(i,:),LineWidth=1);hold on;
% end
% xlabel('Angle of incidence (°)');
% ylabel('Amplitude');
% % ylim([0 0.06]);
% % set(gca,'ytick',[0 0.02 0.04 0.06],'fontsize',Fontsize1);
% grid on
% set(gcf, 'Position', [10 220 500 200]);
time1 = tm.*1e3;
time2 = (thick(1)./Vp2.*1e3).*2 + time1;
time3 = (real(thick(2)./Vp_forward(3,30)).*1e6 + thick(1)./Vp2.*1e3).*2 + time1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.35 scrsz(4)*0.5]); hold on;
wigb(syna_p_r_3(:,1:2:end),scale,ang_inc(1:2:end),tt1);    % 反射率法
% wiggle(syna_p_r_3(:,1:2:end),tt1,ang_inc(1:2:end),1,'d');
yline(time1,'r--',LineWidth=2);
yline(time2,'r--',LineWidth=2);
yline(time3,'r--',LineWidth=2);
ylim([0 350]);
title('(b)','FontName','Times New Roman','Fontsize',Fontsize1);
xlabel('Incident angle (Degree)');
ylabel('Time (ms)');
% set(gcf, 'Position', [50 50 300 500]);
box on;

% =-============ The fourth saturation =====================
% figure,
% for i = 1:length(nd)
%     plot(tt1,syna_p_r_4(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
% end
% ylim([-0.2 0.3]);
% xlim([0 350])
% xlabel('Time (ms)');
% ylabel('amplitude');
% legend('0°','10°','20°','30°');
% grid on
% set(gcf, 'Position', [10 220 500 200]);
% 
% figure;hold on;
% wigb(syna_p_r_4,scale,ang_inc,tt1);    % 反射率法
% ylim([0 350]);
% xlabel('Incident angle (Degree)');
% ylabel('Time');

% =-============ The fifth saturation =====================
% figure,
% for i = 1:length(nd)
%     plot(tt1,syna_p_r_5(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
% end
% ylim([-0.06 0.06]);
% set(gca,'ytick',[-0.04 0 0.04],'fontsize',Fontsize1);
% xlim([0 350])
% xlabel('Time (ms)');
% ylabel('amplitude');
% legend('0°','10°','20°','30°');
% grid on
% set(gcf, 'Position', [10 220 500 200]);

[amp_layer1,index1] = max(abs(syna_p_r_5(1:1000,:)),[],1);
[amp_layer2,index2] = max(abs(syna_p_r_5(1000:2000,:)),[],1);
[amp_layer3,index3] = max(abs(syna_p_r_5(2000:3000,:)),[],1);

% for i = 1:length(amp_layer1)
% syna_p_r_5(:,i) = syna_p_r_5(:,i)./max(abs(syna_p_r_5(:,i)));
% end
% for i = 1:length(index1)
%     amp_layer_1(i) = syna_p_r_5(index1(i),i);
%     amp_layer_2(i) = syna_p_r_5(1000+index2(i),i);
%     amp_layer_3(i) = syna_p_r_5(2000+index3(i),i);
% end
% amp_1ayar = [amp_layer_1 ; amp_layer_2 ; amp_layer_3 ];
% figure;
% for i = 1:3
%     scatter(ang_inc,amp_1ayar(i,:),'color',CC(i,:),LineWidth=1);hold on;
% end
% xlabel('Angle of incidence (°)');
% ylabel('Amplitude');
% % ylim([0 0.06]);
% % set(gca,'ytick',[0 0.02 0.04 0.06],'fontsize',Fontsize1);
% grid on
% set(gcf, 'Position', [10 220 500 200]);
time1 = tm.*1e3;
time2 = (thick(1)./Vp2.*1e3).*2 + time1;
time3 = (real(thick(2)./Vp_forward(5,30)).*1e6 + thick(1)./Vp2.*1e3).*2 + time1;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.35 scrsz(4)*0.5]); hold on;
wigb(syna_p_r_5(:,1:2:end),scale,ang_inc(1:2:end),tt1);     hold on; % 反射率法
% wiggle(syna_p_r_5(:,1:2:end),tt1,ang_inc(1:2:end),1,'d');     hold on;
yline(time1,'r--',LineWidth=2);
yline(time2,'r--',LineWidth=2);
yline(time3,'r--',LineWidth=2);
ylim([0 350]);
title('(c)','FontName','Times New Roman','Fontsize',Fontsize1);
xlabel('Incident angle (Degree)');
ylabel('Time (ms)');
% set(gcf, 'Position', [50 50 300 500]);
box on;

% %  =============  同一道不同饱和度地震记录对比  =================================
% nd1 = 30;
% figure,
% plot(syna_p_r_1(:,nd1),tt1,'b-',LineWidth=1);hold on;
% % plot(syna_p_r_2(:,nd1),tt1,'r-',LineWidth=1);
% plot(syna_p_r_3(:,nd1),tt1,'k',LineWidth=1);
% % plot(syna_p_r_4(:,nd1),tt1,'r-',LineWidth=1);
% plot(syna_p_r_5(:,nd1),tt1,'r',LineWidth=1);
% set(gca,'YDir','reverse');
% ylim([0 450])
% xlabel('amplitude');
% ylabel('Time (ms)');
% legend('Sw = 10%','Sw = 50%','Sw = 90%');
% grid on
% set(gcf, 'Position', [10 220 300 450]);








