clc;    clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----The influence of saturation on the forward simulation results —— Analyze from a fixed Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fontsize = 14;              % For display
AspectRatio = 0.618;
CC1=colormap(jet(7));    CC2=colormap(hot(7));   CC3=colormap(cool(7));   CC4=colormap(spring(7));
red = [1 0 0];  green = [0 1 0];  blue = [0 0 1];   magenta = [1 0 1];   black = [0 0 0];   
CC = [red;blue;black;magenta;green];

load ./model_data./Vp_Sw_perm.mat;
Vp_Sw_perm1 = Vp_crm_1;      % Velocity distribution (varying with saturation and permeability)
[x1,x2,x3] = size(Vp_Sw_perm1);
%%
perm0 = 10^(-15);             %%%%%% Select the permeability        10^(-15) --> 1mD       10^(-12) --> 1D
%%
perm = 10.^(-16:0.01:-12) ;

n_perm = find(ismember(perm,perm0));   
S_1 = 0.01:0.01:0.99;   %  S_1 : saturation of fliud 1

Vp_crm_Sw = reshape(Vp_Sw_perm1(:,n_perm,:),[x1,x3]);

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
% ===============================     Forward simulation    =====================================================  
V_forward = conj(Vp_crm_Sw);
[m,n] = size(V_forward);

Frequency1 = 1:1:200;
values = find(ismember(FrequencyRange,Frequency1));
for i =1:length(Frequency1)
    Vp_forward1(:,i) = V_forward(:,values(i));          
end

%  The velocities and densities of each elastic layer
% Vp1 = 4;
% Vp2 = 4.1;
% Vp4 = 4.1;
% Vs1 = 2.3;
% Vs2 = 2.4;
% Vs4 = 2.4;
% rho1 = 2.3;
% rho2 = 2.4;
% rho4 = 2.4;
Vp1 = 3.9;
Vp2 = 4.1;
Vp4 = 4.1;
Vs1 = 2.1;
Vs2 = 2.3;
Vs4 = 2.3;
rho1 = 2.2;
rho2 = 2.5;
rho4 = 2.5;

rho_s = 2650;
por = 0.12;
rho_f = [rho_water rho_gas];            
rho = rho_s.*(1-por) + por.*S_1.*rho_f(1) + por.*(1 - S_1).*rho_f(2);

Vs_rock = 2200;
for i = 1:m
    Vp_0Hz(i) = min(V_forward(i,:));         
    Vp_forward(i,:) = [Vp_0Hz(i),Vp_forward1(i,:)];    

    Vs_1(i,:) = [Vs1 Vs2 Vs_rock./1000 Vs4];              %  The S-wave velocity of each layer(km/s)

    rho_1(i,:) = [rho1 rho2 rho(i)./1000 rho4];          %  The density of each layer(g/cm^3)
end

Frequency = [0,Frequency1];

thick = [0.2 0.2];          %  The thickness of each layer(km)
nL = length(thick);             %  Number of layers

theta =0;            % Incident angle: 2

for i = 1:m
    for j = 1:length(Frequency)
       Vp_1 = [Vp1,Vp2 , Vp_forward(i,j)./1000 ,Vp4];             %  The P-wave velocity of each layer(km/s)
       [Rm] = ref_multi (Vp_1(1),Vs_1(1),rho_1(1),Vp_1(2:end),Vs_1(i,2:end),rho_1(i,2:end),thick,theta,Frequency(j));
       Rpp_1(i,j) = Rm(1,1);
    end
end

% ----     Frequency domain wavelet   -------------------
f_p = 30;           %  dominant frequency
tm = 0.05;          %  wavelet delay
N = 10000;    
dt = 1./N;
[WW_wavelet]= ricker2(f_p,Frequency,tm);
figure
plot(Frequency,real(WW_wavelet));
xlim([0 100]);

%%%%% ----------------- Seismic synthetics ----============================
for i=1:m
    for j=1:length(Frequency) 
        syna_pp_r(j,i) = (Rpp_1(i,j))*WW_wavelet(j)  ; 
    end
end

for i = 1:m
    [t1,syna_p_r(:,i)] = FFT1(syna_pp_r(:,i),dt,-1);
end

tt1=(0:N-1)*dt*10^3;                 % Time for synthetics, ms

nSw = [0.1 0.5 0.9].*100;         %  Extract several saturation to detect the amplitude differences
nd = find(ismember(S_1.*100,nSw));     
scale = 1;

Fontsize1 = 12;
%  ===================================================================
figure,
for i = 1:length(nd)
    plot(tt1,syna_p_r(:,nd(i)),'color',CC(i,:),LineWidth=1);hold on;
end
% xlim([270 360]);
% ylim([-0.06 0.12]);
xlim([0 350]);
ylim([-0.085 0.085]);
title('(c)','FontName','Times New Roman','Fontsize',Fontsize1);          %  c  or  d
xlabel('Time (ms)');
ylabel('amplitude');
legend('10%','50%','90%');
grid on
set(gca, 'FontSize',Fontsize1);
set(gcf, 'Position', [10 220 500 200]);

figure;
imagesc(S_1.*100,tt1,syna_p_r);
title('(b)','FontName','Times New Roman','Fontsize',Fontsize1);
xlabel('S_W (%)');
ylabel('Time (ms)');
set(gca,'xtick',[10 50 90]);
grid on
ylim([0 350]);
M = seis(2);
colormap(M);
ck=colorbar('horiz');
set(ck,'Position',[0.23 0.2 0.6 0.04]);
clim([-0.1 0.1]);
set(gcf, 'Position', [10 220 300 550]);

