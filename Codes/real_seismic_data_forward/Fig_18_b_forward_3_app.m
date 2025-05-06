clc;
clear;
%%                                        储层分为三层，三层均为常速度

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

%%            井数据分层 ， 速度模型
time(1) = dh./Vp(1);
for i = 2:length(Vp)
    time(i) = dh./Vp(i) + time(i-1);
end
time = 2.*time;  %   双程旅行时
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

%%            合成地震记录
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

d_ang=3;                                      % Angle interval, degree
ang_total = 42;                               % Maximum angle, degree
ang_inc=(3:d_ang:ang_total);                  % Incident angle: 1
theta =(3:d_ang:ang_total)*pi/180;            % Incident angle: 2

Frequency = 0:1:100;

for i = 1:length(theta)
    for j = 1:length(Frequency)
    [Rm] = ref_multi (Vp_model_1(1),Vs_model_1(1),rho_model_1(1),Vp_model_1(2:end),Vs_model_1(2:end),rho_model_1(2:end),thick,theta(i),Frequency(j));
    Rpp(i,j) = Rm(1,1);
    end
end

% ----     频率域子波   -------------------
f_p = 23;          %  主频
tm = 0.125;          %  子波延时
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
wiggle(syna_p_r,tt1,ang_inc,1,'d');   hold on;
yline(tt1(nt1(1)),'b',LineWidth=2);
yline(tt1(nt2(1)),'r',LineWidth=2);
ylim([0 230]);
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Time (ms)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
title('(b)','FontName','Times New Roman','Fontsize',14);
% set(gcf, 'Position', [10 10 600 700]);



% AVA 信息
Amp_Top_3 = syna_p_r(nt1(1),:);
Amp_Base_3 = syna_p_r(nt2(1),:);
% save('./seismic_data./Amp_TOP_3.mat','Amp_Top_3') ;
% save('./seismic_data./Amp_BASE_3.mat','Amp_Base_3') ;

% %  实际数据 与 模拟数据 归一化
% load Amp_base_real.mat;
% load Amp_top_real.mat;
% amp_Top_real = Amp_Top_real./abs(min(Amp_Top_real));
% amp_Base_real = Amp_Base_real./abs(min(Amp_Top_real));
% 
% amp_Top_3 = Amp_Top_3./abs(min(Amp_Top_3));
% amp_Base_3 = Amp_Base_3./abs(min(Amp_Top_3));
% 
% figure
% h1 = scatter(ang_inc,Amp_Top_3.*1000,'blue','filled','o','SizeData',100);    hold on;
% h2 = scatter(ang_inc,Amp_Base_3.*1000,'red','filled','o','SizeData',100);
% h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
% h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
% scatter(ang_inc,amp_Top_3,'blue','filled','o','SizeData',100);    hold on;
% scatter(ang_inc,amp_Base_3,'red','filled','o','SizeData',100);
% scatter(ang_inc,amp_Top_real,'blue','o','SizeData',100,'LineWidth',1.5);
% scatter(ang_inc,amp_Base_real,'red','o','SizeData',100,'LineWidth',1.5);
% % legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
% xlim([5 35]);
% ylim([-2 2]);
% title('(d)','FontName','Times New Roman','Fontsize',14);
% set(gcf, 'Position', [10 10 700 300]);
% set(gca,'fontsize',13);
% grid on; grid minor;  box on;
% legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
% ah=axes('position',get(gca,'position'),'visible','off');
% legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);









