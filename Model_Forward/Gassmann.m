function [Ksat,Gsat] = Gassmann (Kdry,Gdry,Km,fai,Kf)
%利用Gassmann方法加入流体
temp1 = (1 - Kdry./Km).^2;
temp2 = (1 - fai)./Km + fai./Kf - Kdry./(Km.^2);
Ksat = Kdry + temp1./temp2;
Gsat = Gdry;
