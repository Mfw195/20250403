clc;
clear all;
%%            对比各类正演结果

load ./seismic_data./Amp_BASE_real.mat;
load ./seismic_data./Amp_BASE_log.mat;
load ./seismic_data./Amp_BASE_3.mat;            %  储层为3层，常速度
load ./seismic_data./Amp_BASE_3_disp.mat;       %  储层为3层，频散速度

load ./seismic_data./Amp_TOP_real.mat;
load ./seismic_data./Amp_TOP_log.mat;
load ./seismic_data./Amp_TOP_3.mat;            %  储层为3层，常速度
load ./seismic_data./Amp_TOP_3_disp.mat;       %  储层为3层，频散速度

d_ang=3;                                      % Angle interval, degree
ang_total = 42;                               % Maximum angle, degree
ang_inc=(3:d_ang:ang_total);                  % Incident angle

%  归一化
amp_Top_real = Amp_Top_real./max(abs([Amp_Top_real(2:11) Amp_Base_real(2:11)]));
amp_Base_real = Amp_Base_real./max(abs([Amp_Top_real(2:11) Amp_Base_real(2:11)]));

amp_Top_log = Amp_Top_log./max(abs([Amp_Top_log(2:11) Amp_Base_log(2:11)]))+0.2;
amp_Base_log = Amp_Base_log./max(abs([Amp_Top_log(2:11) Amp_Base_log(2:11)]));

amp_Top_3 = Amp_Top_3./max(abs([Amp_Top_3(2:11) Amp_Base_3(2:11)]))+0.1;
amp_Base_3 = Amp_Base_3./max(abs([Amp_Top_3(2:11) Amp_Base_3(2:11)]));

amp_Top_4 = Amp_Top_4./max(abs([Amp_Top_4(2:11) Amp_Base_4(2:11)]))+0.4;
amp_Base_4 = Amp_Base_4./max(abs([Amp_Top_4(2:11) Amp_Base_4(2:11)]));

%      实际数据与测井正演数据
figure
h1 = scatter(ang_inc,amp_Top_log.*1000,'blue','filled','o','SizeData',100);    hold on;
h2 = scatter(ang_inc,amp_Base_log.*1000,'red','filled','o','SizeData',100);
h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
scatter(ang_inc,amp_Top_log,'blue','filled','o','SizeData',100);    hold on;
scatter(ang_inc,amp_Base_log,'red','filled','o','SizeData',100);
scatter(ang_inc,amp_Top_real,'blue','o','SizeData',100,'LineWidth',1.5);    hold on;
scatter(ang_inc,amp_Base_real,'red','o','SizeData',100,'LineWidth',1.5);
% legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlim([5 35]);
ylim([-1.3 2.3]);
title('(c)','FontName','Times New Roman','Fontsize',14);
set(gcf, 'Position', [10 100 700 300]);
set(gca,'fontsize',13);
grid on; grid minor;  box on;
legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);

%      实际数据与简化速度模型正演数据
figure
h1 = scatter(ang_inc,amp_Top_3.*1000,'blue','filled','o','SizeData',100);    hold on;
h2 = scatter(ang_inc,amp_Base_3.*1000,'red','filled','o','SizeData',100);
h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
scatter(ang_inc,amp_Top_3,'blue','filled','o','SizeData',100);    hold on;
scatter(ang_inc,amp_Base_3,'red','filled','o','SizeData',100);
scatter(ang_inc,amp_Top_real,'blue','o','SizeData',100,'LineWidth',1.5);
scatter(ang_inc,amp_Base_real,'red','o','SizeData',100,'LineWidth',1.5);
% legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlim([5 35]);
ylim([-1.3 2.3]);
title('(d)','FontName','Times New Roman','Fontsize',14);
set(gcf, 'Position', [10 10 700 300]);
set(gca,'fontsize',13);
grid on; grid minor;  box on;
legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);

%      实际数据与频散速度模型正演数据
figure
h1 = scatter(ang_inc,amp_Top_4.*1000,'blue','filled','o','SizeData',100);    hold on;
h2 = scatter(ang_inc,amp_Base_4.*1000,'red','filled','o','SizeData',100);
h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
scatter(ang_inc,amp_Top_4,'blue','filled','o','SizeData',100);    hold on;
scatter(ang_inc,amp_Base_4,'red','filled','o','SizeData',100);
scatter(ang_inc,amp_Top_real,'blue','o','SizeData',100,'LineWidth',1.5);
scatter(ang_inc,amp_Base_real,'red','o','SizeData',100,'LineWidth',1.5);
% legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlim([5 35]);
ylim([-1.3 2.3]);
title('(d)','FontName','Times New Roman','Fontsize',14);
set(gcf, 'Position', [10 10 700 300]);
set(gca,'fontsize',13);
grid on; grid minor;  box on;
legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);

real_matrix = [amp_Top_real(:,2:11)',amp_Base_real(:,2:11)'];
log_matrix = [amp_Top_log(:,2:11)',amp_Base_log(:,2:11)'];
Elasticity_matrix = [amp_Top_3(:,2:11)',amp_Base_3(:,2:11)'];
Viscoelasticity_matrix = [amp_Top_4(:,2:11)',amp_Base_4(:,2:11)'];

corr_1 = corrcoef(real_matrix,log_matrix);
corr_2 = corrcoef(real_matrix,Elasticity_matrix);
corr_3 = corrcoef(real_matrix,Viscoelasticity_matrix);             % 相关系数

norm_1 = norm(real_matrix-log_matrix)/norm(real_matrix);
norm_2 = norm(real_matrix-Elasticity_matrix)/norm(real_matrix);
norm_3 = norm(real_matrix-Viscoelasticity_matrix)/norm(real_matrix);  % 相对误差

%%
%      弹性速度模型与频散速度模型正演数据
figure
h1 = scatter(ang_inc,amp_Top_4.*1000,'blue','filled','o','SizeData',100);    hold on;
h2 = scatter(ang_inc,amp_Base_4.*1000,'red','filled','o','SizeData',100);
h3 = scatter(ang_inc,amp_Top_real.*1000,'k','o','SizeData',100,'LineWidth',1.5);
h4 = scatter(ang_inc,amp_Base_real.*1000,'k','filled','o','SizeData',100,'LineWidth',1.5);         %  for legend
scatter(ang_inc,amp_Top_4,'blue','filled','o','SizeData',100);    hold on;
scatter(ang_inc,amp_Base_4,'red','filled','o','SizeData',100);
scatter(ang_inc,amp_Top_3,'blue','o','SizeData',100,'LineWidth',1.5);
scatter(ang_inc,amp_Base_3,'red','o','SizeData',100,'LineWidth',1.5);
% legend('Top of the reservoir','Base of the reservoir','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
xlim([5 35]);
ylim([-2 2]);
title('(d)','FontName','Times New Roman','Fontsize',14);
set(gcf, 'Position', [10 10 700 300]);
set(gca,'fontsize',13);
grid on; grid minor;  box on;
legend([h1 h2],'Top of the reservoir','Base of the reservoir','Fontsize',13);
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,[h3 h4],'Real seismic data','Synthetic data','Fontsize',13);