function [ff1]= ricker2(fm,f,t0)
% �׿��Ӳ�����Ƶ����ı��ʽ
%fm:��Ƶ��
%f:Ƶ��
ff1=2.*(f).^2.*exp(-(f/fm).^2).*exp(-1j.*2.*pi.*f.*t0)/(fm^3.*sqrt(pi));
end

