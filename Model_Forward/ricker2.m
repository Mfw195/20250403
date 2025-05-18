function [ff1]= ricker2(fm,f,t0)
% 雷克子波关于频率域的表达式
%fm:主频率
%f:频率
ff1=2.*(f).^2.*exp(-(f/fm).^2).*exp(-1j.*2.*pi.*f.*t0)/(fm^3.*sqrt(pi));
end

