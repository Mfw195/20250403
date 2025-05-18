clc;
clear all;
%%                          Wellbore bypass seismic data are synthesized from logging data

% load ./log_data/depth_log.mat;
% load ./log_data/rho_log.mat;
% load ./log_data/Vp_log.mat;
% load ./log_data/Vs_log.mat;

load ./log_data_new/depth_log.mat;
load ./log_data_new/rho_log.mat;
load ./log_data_new/Vp_log.mat;
load ./log_data_new/Vs_log.mat;

dh = depth(2)-depth(1);
thick = dh.*ones(1,length(Vp)-2).*0.001;


d_ang=3;                                      % Angle interval, degree
ang_total = 42;                               % Maximum angle, degree
ang_inc=(3:d_ang:ang_total);                  % Incident angle: 1
theta =(3:d_ang:ang_total)*pi/180;            % Incident angle: 2

Frequency = 0:1:100;

for i = 1:length(theta)
    for j = 1:length(Frequency)
    [Rm] = ref_multi (Vp(1),Vs(1),rho(1),Vp(2:end),Vs(2:end),rho(2:end),thick,theta(i),Frequency(j));
    Rpp(i,j) = Rm(1,1);
    end
end

% ----     Frequency domain wavelet   -------------------
f_p = 23;          %  dominant frequency
tm = 0;          %  Wavelet delay
N = 10000;    
dt = 1./N;
[WW_wavelet]= ricker2(f_p,Frequency,tm);
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

scale = 1;
figure;
wiggle(syna_p_r,tt1,ang_inc,1,'d');
yline(tt1(nt1),'b',LineWidth=1);
yline(tt1(nt2),'r',LineWidth=1);
ylim([0 230]);
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Time (ms)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
title('(a)','FontName','Times New Roman','Fontsize',14);
% set(gcf, 'Position', [10 10 600 700]);

figure;
wigb(syna_p_r(1:2301,3:4),scale,ang_inc(3:4),tt1(1:2301));     hold on;
yline(tt1(nt1),'b',LineWidth=2);
yline(tt1(nt2),'r',LineWidth=2);
ylim([0 230]);
xlabel('Incident angle (Degree)');
ylabel('Time');
set(gcf, 'Position', [10 10 300 500]);

Amp_Top_log = syna_p_r(nt1(1),:);
Amp_Base_log = syna_p_r(nt2(1),:);
% save('./seismic_data./Amp_TOP_log.mat','Amp_Top_log') ;
% save('./seismic_data./Amp_BASE_log.mat','Amp_Base_log') ;

% load Amp_base_real.mat;
% load Amp_top_real.mat;
% 


% amp_Top_real = Amp_Top_real./max(abs([Amp_Top_real(2:11)]));
% amp_Base_real = Amp_Base_real./abs(max([Amp_Base_real(2:11) Amp_Top_real(2:11)]));
% 
% amp_Top_log = Amp_Top_log./abs(max([Amp_Top_log(2:11) Amp_Base_log(2:11)]));
% amp_Base_log = Amp_Base_log./abs(max([Amp_Top_log(2:11) Amp_Base_log(2:11)]));
% 
% figure
% h1 = scatter(ang_inc,Amp_Top_log.*1000,'blue','filled','o','SizeData',100);    hold on;
% h2 = scatter(ang_inc,Amp_Base_log.*1000,'red','filled','o','SizeData',100);
% h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
% h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
% scatter(ang_inc,amp_Top_log,'blue','filled','o','SizeData',100);    hold on;
% scatter(ang_inc,amp_Base_log,'red','filled','o','SizeData',100);
% scatter(ang_inc,amp_Top_real,'blue','o','SizeData',100,'LineWidth',1.5);    hold on;
% scatter(ang_inc,amp_Base_real,'red','o','SizeData',100,'LineWidth',1.5);
% % legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% xlim([5 35]);
% ylim([-2 2]);
% title('(c)','FontName','Times New Roman','Fontsize',14);
% set(gcf, 'Position', [10 100 700 300]);
% set(gca,'fontsize',13);
% grid on; grid minor;  box on;
% legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
% ah=axes('position',get(gca,'position'),'visible','off');
% legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);

