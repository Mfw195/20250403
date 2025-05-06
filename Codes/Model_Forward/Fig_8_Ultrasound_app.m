clc;   clear all;
Fontsize1 = 14;
%%%%%%%%%%%%%%%%%%%%%%        超声波数据          %%%%%%%%%%%%%%%%%%
Press = [0 1 3 5 8 10 15 20 25 30 35 40 45 50];

load ./Experimental_data./Ultrasound_data./VP_dry;
load ./Experimental_data./Ultrasound_data./VS_dry;

load ./Experimental_data./Ultrasound_data./VP_water;
load ./Experimental_data./Ultrasound_data./Vs_water;

load ./Experimental_data./Ultrasound_data./VP_glycerol;
load ./Experimental_data./Ultrasound_data./VS_glycerol;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(Press, Vp_dry./1000,'-kd','MarkerFaceColor','k','MarkerSize',7,LineWidth=1.5);  hold on;
plot(Press, Vp_water./1000,'-ro','MarkerFaceColor','r','MarkerSize',7,LineWidth=1.5);
plot(Press, Vp_ganyou./1000,'-b^','MarkerFaceColor','b','MarkerSize',7,LineWidth=1.5);
xlabel('Effective pressure {\itP_{eff}} (MPa)','Fontsize',Fontsize1);
ylabel('{\itV}_P (km/s)','Fontsize',Fontsize1); 
set(gca,'xtick',[0 10 20 30 40 50],'Fontsize',Fontsize1);
set(gca,'ytick',[3 3.5 4 4.5 5],'Fontsize',Fontsize1);
legend('Dry','Brine','Glycerol');
ylim([2.6 5.1]);
xlim([-2 52]);
title('(a)','FontName','Times New Roman','Fontsize',Fontsize1);
grid on;   box on;  grid minor;


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(Press, Vs_dry./1000,'-kd','MarkerFaceColor','k','MarkerSize',7,LineWidth=1.5);  hold on;
plot(Press, Vs_water./1000,'-ro','MarkerFaceColor','r','MarkerSize',7,LineWidth=1.5);
plot(Press, Vs_ganyou./1000,'-b^','MarkerFaceColor','b','MarkerSize',7,LineWidth=1.5);
xlabel('Effective pressure {\itP_{eff}} (MPa)','Fontsize',Fontsize1);
ylabel('{\itV}_S (km/s)','Fontsize',Fontsize1); 
set(gca,'xtick',[0 10 20 30 40 50],'Fontsize',Fontsize1);
set(gca,'ytick',[2 2.25 2.5 2.75 3],'Fontsize',Fontsize1);
legend('Dry','Brine','Glycerol');
ylim([2 3]);
xlim([-2 52]);
title('(b)','FontName','Times New Roman','Fontsize',Fontsize1);
grid on;   box on;  grid minor;















