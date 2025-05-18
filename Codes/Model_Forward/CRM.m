function [Hh,Hw,H_sc,V,Q] = CRM(Kdry1,Kdry2,Gdry1,Gdry2,Ks,fai,S1,S2,Kf1,Kf2,vis_f1,vis_f2,rho,omiga,perm,b,xuanze)
%   三维随机斑块饱和 
%   xuanze = 1  -- > 指数型相关函数
%   xuanze = 2  -- > 高斯型相关函数
%%%%%%%%%%%%%%%%%%%%%%%       求Gassmann Hill和Gassmann Wood的纵波模量%%%%%%%%%%%%
%  S1,kf1空气，            S2,kf2水

Gsat1 = Gdry1;   Gsat2 = Gdry2;
Gsat = S1.*Gsat1 + Gsat2.*S2;

alpha1 = 1-Kdry1./Ks;            alpha2 = 1-Kdry2./Ks;          %Biot系数

Kdry=S1.*Kdry1+S2.*Kdry2;          Gdry=S1.*Gdry1+S2.*Gdry2;    % 湿骨架体积模量剪切模量                     
alpha=1-Kdry./Ks;        

M1 = 1./((alpha1-fai)./Ks + fai./Kf1);    M2 = 1./((alpha2-fai)./Ks + fai./Kf2);                         
% 孔隙空间模量                                              
Ldry = Kdry + 4*Gdry/3;         

%干骨架纵波模量
Ksat1 = Kdry1 + alpha1.^2.*M1;
Ksat2 = Kdry2 + alpha2.^2.*M2;
Lsat1 = Ksat1 + 4*Gdry1/3;
Lsat2 = Ksat2 + 4*Gdry2/3;

Hh = 1./(S1./Lsat1 + S2./Lsat2);             % equ(3) Gassmann Hill

Kf = 1./(S1./Kf1+S2./Kf2);                                 % equ(1) S1空气，S2水
M = 1./((alpha-fai)./Ks+fai./Kf);                          % equ(7a,b)下面那一段
Ksat = Kdry+alpha.^(2).*M;
Hw = Ksat + Gdry*4/3;                                  % equ(1) Gassmann Wood

M0=M1.*S1+M2.*S2;                                        % equ(9)
H0=Ldry+alpha.^(2).*M0;                                   % equ(7a,b)




%%%%%%%%%%%%%%%%%%%%%%%       扰动量      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cha = M2-M1;                     
x_1 = -cha/2;                                   %%%扰动量
x_2 = cha/2;                          
X1 = M1 + cha/2;                 
epsion_m1 = x_1./X1;
epsion_m2 = x_2./X1;
delta_mm = epsion_m1.^(2)*S1+epsion_m2.^(2)*S2;   %  Muller and Gurevich,2005a,equ(19)(影响衰减初始值大小，正比)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta2 = alpha.^(2).*M0.*delta_mm.^(2)./2./H0;                                % equ(8a,b)
delta1 = Ldry.*delta2./H0;
N = M0.*Ldry./H0;                                                          % equ(6b)
vis = vis_f1.*S1 + vis_f2.*S2; %混和流体粘度                                  % equ(9)
k_ = sqrt(1i*omiga*vis./perm./N);    %慢波数                             % equ(6a,b)
%k_jia = omiga*sqrt(rho/H0);         %快波数

H_low = H0.*(1-delta2).^2;                                                % equ(13)
H_high = H0.*(1-delta2+delta1).^2;                                        % equ(14)
if xuanze == 1
    H_eff = H0.*(1-delta2-(delta1.*k_.^(2)*b^(2))./(1i*k_.*b-1).^2).^2;       % equ(16)
else
    y = -1i*k_.*b./2;
    erfcy = double(feval(symengine,'erfc', y));
    H_eff = H0.*(1-delta2 + 2*delta1.*y.^2.*(1 - y.*sqrt(pi).*exp(y.^2).*erfcy)).^2;       % equ(17)
end
H_sc = Hw.*(1+(Hh-Hw)./(H_high-H_low).*(H_eff-H_low)./Hw);                  % equ(15)
V = sqrt((H_sc)/rho);                                               % equ(11)
Q = abs(imag(H_sc))./real(H_sc);
% Vs = sqrt(Gsat/rho);
% Qs = abs(imag(Gsat))./real(Gsat);
% save 'V_crm.mat' V;
% save 'Q_crm.mat' Q;