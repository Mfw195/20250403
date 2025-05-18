clc;
clear;
%%     The reservoir is divided into three layers, and all three layers are of dispersion velocity

load ./log_data_new/depth_log.mat;
load ./log_data_new/rho_log.mat;
load ./log_data_new/Vp_log.mat;
load ./log_data_new/Vs_log.mat;

dh = depth(2)-depth(1);
thick = dh.*ones(1,length(Vp)-2).*0.001;

nt1 = find(fix(depth) == 3308);
nt2 = find(fix(depth) == 3339);
nt3 = find(fix(depth) == 3356);
nt4 = find(fix(depth) == 3374);
nt5 = find(fix(depth) == 3530);

%%           Well data stratification is used as the velocity model
time(1) = dh./Vp(1);
for i = 2:length(Vp)
    time(i) = dh./Vp(i) + time(i-1);
end
time = 2.*time; 
time_up = time(1:nt1(1));
time_rev_1 = time(nt1(2):nt2(1));
time_rev_2 = time(nt2(2):nt3(1));
time_rev_3 = time(nt3(2):nt4(1));
time_down = time(nt4(2):end);

Vp_up = mean(Vp(1:nt1(1))).*ones(1,length(time_up));
Vp_rev_1 = mean(Vp(nt1(2):nt2(1))).*ones(1,length(time_rev_1));
Vp_rev_2 = mean(Vp(nt2(2):nt3(1))).*ones(1,length(time_rev_2));
Vp_rev_3 = mean(Vp(nt3(2):nt4(1))).*ones(1,length(time_rev_3));
Vp_down = mean(Vp(nt4(2):nt5(end))).*ones(1,length(time_down));

rho_up = mean(rho(1:nt1(1))).*ones(1,length(time_up));
rho_rev_1 = mean(rho(nt1(2):nt2(1))).*ones(1,length(time_rev_1));
rho_rev_2 = mean(rho(nt2(2):nt3(1))).*ones(1,length(time_rev_2));
rho_rev_3 = mean(rho(nt3(2):nt4(1))).*ones(1,length(time_rev_3));
rho_down = mean(rho(nt4(2):nt5(end))).*ones(1,length(time_down));

rho_model = [rho_up rho_rev_1 rho_rev_2 rho_rev_3 rho_down];
Vp_model = [Vp_up Vp_rev_1 Vp_rev_2 Vp_rev_3 Vp_down];

thick_up = Vp_up(1).*time_up(end)./1e3;
thick_rev_1 = Vp_rev_1(1).*(time_rev_1(end) - time_up(end))./1e3;
thick_rev_2 = Vp_rev_2(1).*(time_rev_2(end) - time_rev_1(end))./1e3;
thick_rev_3 = Vp_rev_3(1).*(time_rev_3(end) - time_rev_2(end))./1e3;
thick_down = Vp_down(1).*(time_down(end) - time_rev_3(end))./1e3;

rho_model_1 = [rho_up(1) rho_rev_1(1) rho_rev_2(1) rho_rev_3(1) rho_down(1)];
Vp_model_1 = [Vp_up(1) Vp_rev_1(1) Vp_rev_2(1) Vp_rev_3(1) Vp_down(1)];
rand_1 = rand(length(Vp_model_1),1)';
Vs_model_1 = Vp_model_1.*0.583 - 0.078;

thick = [ thick_rev_1 thick_rev_2 thick_rev_3]./2;
%%                   Synthetic seismic record
% FrequencyRange=[10.^(-4:0.01:-0.01) 1:1:200 10.^(2.4:0.01:8)];
% Frequency1 = 1:1:101;
% values = find(ismember(FrequencyRange,Frequency1));
% load Vp_perm.mat;
% load Vp_Sw.mat;
% 
% Vp_perm = Vp_crm_perm(:,values)./1e3;
% Vp_Sw = Vp_crm_Sw(:,values)./1e3;

d_ang=3;                                      % Angle interval, degree
ang_total = 42;                               % Maximum angle, degree
ang_inc=(3:d_ang:ang_total);                  % Incident angle: 1
theta =(3:d_ang:ang_total)*pi/180;            % Incident angle: 2
Frequency = 0:1:100;
 
% Vp_pinsan_1 = Vp_perm(5,:);
% Vp_pinsan_2 = Vp_Sw(5,:);
% Vp_pinsan_3 = Vp_perm(3,:);
% save('./seismic_data./Vp_1.mat','Vp_pinsan_1') ; 
% save('./seismic_data./Vp_2.mat','Vp_pinsan_2') ; 
% save('./seismic_data./Vp_3.mat','Vp_pinsan_3') ; 

load seismic_data\Vp_1.mat;
load seismic_data\Vp_2.mat;
load seismic_data\Vp_3.mat;

scrsz = get(0,'ScreenSize');   
figure('Position',[400 400 scrsz(4)*0.55 scrsz(4)*0.32]);
plot(Frequency,Vp_pinsan_1,'r-','linewidth',1.5);   hold on;
plot(Frequency,Vp_pinsan_2,'b-','linewidth',1.5);
plot(Frequency,Vp_pinsan_3,'k-','linewidth',1.5);
legend('Thin layar 1','Thin layar 2','Thin layar 3','Fontsize',10);
title('(a)','FontName','Times New Roman','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('P-wave velocity (km/s)','Fontsize',14);
set(gca,'xtick',[0 20 40 60 80 100],'Fontsize',14);
set(gca,'ytick',[3.7 3.8 3.9 4],'Fontsize',14);
grid on; grid minor; box on;

Qp_pinsan_1 = 2.*abs(imag(Vp_pinsan_1)./real(Vp_pinsan_1));
Qp_pinsan_2 = 2.*abs(imag(Vp_pinsan_2)./real(Vp_pinsan_2));
Qp_pinsan_3 = 2.*abs(imag(Vp_pinsan_3)./real(Vp_pinsan_3));
scrsz = get(0,'ScreenSize');   
figure('Position',[400 400 scrsz(4)*0.55 scrsz(4)*0.32]);
plot(Frequency,Qp_pinsan_1,'r-','linewidth',1.5);   hold on;
plot(Frequency,Qp_pinsan_2,'b-','linewidth',1.5);
plot(Frequency,Qp_pinsan_3,'k-','linewidth',1.5);
legend('Thin layar 1','Thin layar 2','Thin layar 3','Fontsize',10);
title('(b)','FontName','Times New Roman','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',14);
ylabel('1/Q_p','Fontsize',14);
ylim([0 0.06]);
set(gca,'xtick',[0 20 40 60 80 100],'Fontsize',14);
set(gca,'ytick',[0  0.02  0.04 0.06 0.08 0.1],'Fontsize',14);
grid on; grid minor; box on;

for i = 1:length(theta)
    for j = 1:length(Frequency)
        Vp_model_1 = [Vp_up(1) Vp_pinsan_1(j) Vp_pinsan_2(j) Vp_pinsan_3(j) Vp_down(1)];
        rand_1 = rand(length(Vp_model_1),1)';
        Vs_model_1 = Vp_model_1.*0.57 - 0.078;
    [Rm] = ref_multi (Vp_model_1(1),Vs_model_1(1),rho_model_1(1),Vp_model_1(2:end),Vs_model_1(2:end),rho_model_1(2:end),thick,theta(i),Frequency(j));
    Rpp(i,j) = Rm(1,1);
    end
end

% ----     Frequency domain wavelet  -------------------
f_p = 23;          %  dominant frequency
tm = 0.125;          %  Wavelet delay
N = 10000;    
dt = 1./N;
[WW_wavelet]= ricker2(f_p,Frequency,tm);
figure
plot(Frequency,real(WW_wavelet));
xlim([0 100]);

%%%%% ----------------- Seismic synthetics ----============================
for i=1:length(theta)
    for j=1:length(Frequency) 
        syna_pp_r(j,i) = (Rpp(i,j))*WW_wavelet(j)  ;
    end
end

for i = 1:length(theta)
    [t1,syna_p_r(:,i)] = FFT1(syna_pp_r(:,i),dt,-1);
end

tt1=(0:N-1)*dt*10^3;                 % Time for synthetics, ms
nt1 = find(fix(tt1) == 125);
nt2 = find(fix(tt1) == 161);

figure;
wiggle(syna_p_r,tt1,ang_inc,1,'d');  
yline(tt1(nt1(1)),'b',LineWidth=2);
yline(tt1(nt2(1)),'r',LineWidth=2);
% xlim([5 35]);
title('(a)','FontName','Times New Roman','Fontsize',14);
ylim([0 230]);
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Time (ms)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% set(gcf, 'Position', [10 10 600 700]);

Amp_Top_4 = syna_p_r(nt1(1),:);
Amp_Base_4 = syna_p_r(nt2(1),:);
% save('./seismic_data./Amp_TOP_3_disp.mat','Amp_Top_4') ;
% save('./seismic_data./Amp_BASE_3_disp.mat','Amp_Base_4') ;


% figure
% scatter(ang_inc,Amp_Top_4,'blue','filled','o','LineWidth',2);    hold on;
% scatter(ang_inc,Amp_Base_4,'red','filled','o','LineWidth',2);
% legend('Top of the reservoir','Base of the reservoir');
% xlabel('Incident angle (Degree)');
% ylabel('Amplitude');
% xlim([5 35]);
% set(gcf, 'Position', [10 10 700 300]);
% set(gca,'fontsize',13);
% grid on;   box on;






























