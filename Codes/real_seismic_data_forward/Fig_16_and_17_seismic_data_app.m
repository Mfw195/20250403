clc;
clear all;
%%     叠后地震剖面
Data= altreadsegy('arbitrary_line-5-7_0_12.sgy');
dao = 250;
[m,n] = size(Data);

% t = linspace(1400,2300,m);
% tt = t(306:511).*1.122-400;
% save('./seismic_data./tt_seismic.mat','tt') ;         %   时间已保存

load seismic_data\tt_seismic.mat;
% data = Data(306:511,:);
data = Data(311:515,:);

figure
imagesc([1 n],tt,data(:,200:300)); hold on;
plot([250,250],get(gca, 'YLim'), '-b', 'LineWidth', 2);
xlabel('CMP','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Time (ms)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
h =  colorbar;
M = seis(2);
colormap(M);
title('(a)','FontName','Times New Roman','Fontsize',14);
set(get(h,'Title'),'string','Amplitude','FontName','Times New Roman','Fontsize',11);
clim;
txt = '\leftarrow well 119-13';
text(250,1550,txt,'FontSize',14,'FontWeight','bold');

%%              实际叠前地震数据
load ./pre_gather.mat;

ang_inc = 3:3:42;
% t = linspace(1400,2300,901);
% tt = t(306:511).*1.122-200;
% pre_data = pre_gather(306:511,:);
pre_data = pre_gather(311:516,:);

t2 = tt - tt(1);
ntt1 = find(fix(t2) == 125);
ntt2 = find(fix(t2) == 161);

figure
wiggle(pre_data,tt,ang_inc,1,'d');
yline(tt(ntt1),'b',LineWidth=1.5);
yline(tt(ntt2),'r',LineWidth=1.5);
set(gca,'xtick',[0 10 20 30 40]);
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Time (ms)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
title('(b)','FontName','Times New Roman','Fontsize',14);

Amp_Top_real = pre_data(ntt1(1),:);
Amp_Base_real = pre_data(ntt2(1),:);
save('./seismic_data./Amp_TOP_real.mat','Amp_Top_real') ;
save('./seismic_data./Amp_BASE_real.mat','Amp_Base_real') ;

figure
scatter(ang_inc,pre_data(ntt1,:),'blue','filled','o','LineWidth',2);    hold on;
scatter(ang_inc,pre_data(ntt2,:),'red','filled','o','LineWidth',2);
legend('Top of the reservoir','Base of the reservoir');
xlim([5 35]);
xlabel('Incident angle (Degree)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Amplitude','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gcf, 'Position', [10 10 700 300]);
set(gca,'fontsize',13);
grid on; grid minor;   box on;
%%     测井数据

%  原始井数据
% load ./log_data/depth_log.mat;
% load ./log_data/rho_log.mat;
% load ./log_data/Vp_log.mat;
% load ./log_data/Vs_log.mat;

%  修改井数据
load ./log_data_new/depth_log.mat;
load ./log_data_new/rho_log.mat;
load ./log_data_new/Vp_log.mat;
load ./log_data_new/Vs_log.mat;

dh = depth(2) - depth(1);
time(1) = dh./Vp(1);
for i = 2:length(Vp)
    time(i) = dh./Vp(i) + time(i-1);
end
time = 2.*time;  %   双程旅行时

% nt1 = find(fix(time) == 131);
% nt2 = find(fix(time) == 148);
% nt3 = find(fix(time) == 156);
% nt4 = find(fix(time) == 165);
nt5 = find(fix(time) == fix(t2(end)));


%  测井数据 
figure
plot(rho(1:8:end),depth(1:8:end),'b',LineWidth=1.5);  hold on;
xlabel('Density (g/cm^3)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
title('(a)','FontName','Times New Roman','Fontsize',14);
set(gca,'YDir','reverse');
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);
%  隐藏Y轴
% set(gca,'yticklabel','');

figure
plot(Vp(1:8:end),depth(1:8:end),'b',LineWidth=1.5);  hold on;
xlabel('{\itV}_P (km/s)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gca,'YDir','reverse');
title('(b)','FontName','Times New Roman','Fontsize',14);
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);
%  隐藏Y轴
% set(gca,'yticklabel','');

figure
plot(Vs(1:8:end),depth(1:8:end),'b',LineWidth=1.5);  hold on;
xlabel('{\itV}_S (km/s)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gca,'YDir','reverse');
title('(c)','FontName','Times New Roman','Fontsize',14);
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);
%  隐藏Y轴
% set(gca,'yticklabel','');


%%      井数据分层

nt1 = find(fix(depth) == 3308);
nt2 = find(fix(depth) == 3339);
nt3 = find(fix(depth) == 3356);
nt4 = find(fix(depth) == 3374);
nt5 = find(fix(depth) == 3530);

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

Al_up = Vp_up.*rho_up;
Al_rev_1 = Vp_rev_1.*rho_rev_1;
Al_rev_2 = Vp_rev_2.*rho_rev_2;
Al_rev_3 = Vp_rev_3.*rho_rev_3;
Al_down = Vp_down.*rho_down;

rho_model = [rho_up rho_rev_1 rho_rev_2 rho_rev_3 rho_down];
Vp_model = [Vp_up Vp_rev_1 Vp_rev_2 Vp_rev_3 Vp_down];
Al_model = [Al_up Al_rev_1 Al_rev_2 Al_rev_3 Al_down];
Vs_model = Vp_model.*0.583 - 0.078;

figure
plot(rho_model,depth,'b',LineWidth=1.5);  hold on;
xlabel('Density (g/cm^3)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gca,'YDir','reverse');
title('(d)','FontName','Times New Roman','Fontsize',14);
set(gca,'YDir','reverse');
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);

figure
plot(Vp_model,depth,'b',LineWidth=1.5);  hold on;
xlabel('{\itV}_P (km/s)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gca,'YDir','reverse');
title('(e)','FontName','Times New Roman','Fontsize',14);
set(gca,'YDir','reverse');
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);

figure
plot(Vs_model,depth,'b',LineWidth=1.5);  hold on;
xlabel('{\itV}_S (km/s)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
ylabel('Depth (m)','FontName','Times New Roman','Fontsize',15,'FontWeight','bold');
set(gca,'YDir','reverse');
title('(f)','FontName','Times New Roman','Fontsize',14);
set(gca,'YDir','reverse');
ylim([depth(1) depth(nt5(10))]);
set(gcf, 'Position', [10 20 200 700]);







