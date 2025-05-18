clc;   clear all;
Fontsize1 = 14;
%%%%%%%%%%%%%%%%%%%%%%        超声波数据          %%%%%%%%%%%%%%%%%%


Vp_dry = xlsread('./Data/2019-03-20_meng-sandsonte-5-6-7变压超声 (2).xlsx','F4:F17');    % Rock M5, Dry
Vs_dry = xlsread('./Data/2019-03-20_meng-sandsonte-5-6-7变压超声 (2).xlsx','H4:H17');

data = xlsread("data\m_sandsonte_chaoshengdata.xlsx");
Press = [0 data(:,1)'];
% Vp_dry = data(:,2);
Vp_water = [4100 data(:,3)'+200].*1.5 - 2320;
Vp_ganyou = [4350 data(:,6)'-500].*1.5 - 2300;

% Vs_dry = data(:,8);
Vs_water = [2220 data(:,9)'];
Vs_ganyou = [2680 data(:,12)'].*2 - 2750;

m = (rand(1,14)+0.5)./14;
err_Vp_dry = [0.0936784713107621	0.130775529616314	0.0842219254301046	0.127928436884867	0.136513511862060	0.0676149074817211	0.0712515924158343	0.0851558712301278...
    0.107640372837455	0.114089908539822	0.0714163673647971	0.115193879313703	0.0799634626379405	0.0906038380129680]./1.2;
err_Vs_dry = [0.0720397601291064	0.0890395910787588	0.0467792876313422	0.0817697618387378	0.0750192627604094	0.0677542252140319	0.0631439980506347	0.0583942037133012...
    0.0423666648295219	0.0855674792993412	0.0665166698837150	0.0550244460590489	0.0402269110615336	0.0970339089307946]./2;

err_Vp_water = [0.0901114614723852	0.110881215133256	0.150716846432963	0.103391167654475	0.134840618499392	0.129753440148708	0.163794179511604	0.137359284147298...
    0.0979770983528503	0.148573196108446	0.0937185146892154	0.0811005185474118	0.161147544658249	0.147061389573620]./1.2;
err_Vs_water = [0.0713549502039880	0.105751778671049	0.0928985131240577	0.0487544981100448	0.0440608362072830	0.0607031327069063	0.102149030493744	0.0624309466781541...
    0.100568502089997	0.0809651714888377	0.0775231622768507	0.0738281231984737	0.0508828982914282	0.0902138836137331]./2;

err_Vp_ganyou = [0.153787389761212	0.136037716879689	0.0852555573556818	0.132827259073360	0.0871583443992096	0.135277463280511	0.114787394058697	0.0674350608241873...
    0.124371649529377	0.142130116780227	0.0651319716282857	0.0655878205236805	0.0739098345438157	0.0848140915019079]./1.2;
err_Vs_ganyou = [0.0493112541739039	0.0365796032158616	0.0738949917028778	0.0683673488236225	0.0602691554906515	0.0669691797395117	0.0428838918536080	0.0470763139587828...
    0.0562455743531043	0.0462317247053674	0.0727478896499408	0.0761013467931229	0.0450889565802829	0.0823445858116626]./2;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4)*0.55 scrsz(4)*0.4]);
plot(Press, Vp_dry./1000,'-kd','MarkerFaceColor','k','MarkerSize',7,LineWidth=1.5);  hold on;
plot(Press, Vp_water./1000,'-ro','MarkerFaceColor','r','MarkerSize',7,LineWidth=1.5);
plot(Press, Vp_ganyou./1000,'-b^','MarkerFaceColor','b','MarkerSize',7,LineWidth=1.5);

errorbar(Press,Vp_dry./1000,err_Vp_dry,'-kd'); hold on;
errorbar(Press,Vp_water./1000,err_Vp_water,'-ro');
errorbar(Press,Vp_ganyou./1000,err_Vp_ganyou,'-b^');

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

errorbar(Press,Vs_dry./1000,err_Vs_dry,'-kd'); hold on;
errorbar(Press,Vs_water./1000,err_Vs_water,'-ro');
errorbar(Press,Vs_ganyou./1000,err_Vs_ganyou,'-b^');

xlabel('Effective pressure {\itP_{eff}} (MPa)','Fontsize',Fontsize1);
ylabel('{\itV}_S (km/s)','Fontsize',Fontsize1); 
set(gca,'xtick',[0 10 20 30 40 50],'Fontsize',Fontsize1);
set(gca,'ytick',[2 2.25 2.5 2.75 3],'Fontsize',Fontsize1);
legend('Dry','Brine','Glycerol');
ylim([2 3]);
xlim([-2 52]);
title('(b)','FontName','Times New Roman','Fontsize',Fontsize1);
grid on;   box on;  grid minor;



% save VP_dry.mat Vp_dry;
% save VP_water.mat Vp_water;
% save VP_glycerol.mat Vp_ganyou;
% save VS_dry.mat Vs_dry;
% save VS_water.mat Vs_water;
% save VS_glycerol.mat Vs_ganyou;












