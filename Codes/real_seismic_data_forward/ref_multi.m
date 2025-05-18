function [Rm] = ref_multi (Vp0,Vs0,rho0,Vp,Vs,rho,d,seta,f)
%  Vp0,Vs0,rho0 :   第0层速度密度
%  Vp,Vs,rho    :   第1-n 层速度密度

n = length(Vp);
R1 = zeros(2,2,n-1);
R2 = zeros(2,2,n-1);
Rd = zeros(2,2,n);
Ru = zeros(2,2,n-1);
Td = zeros(2,2,n-1);
Tu = zeros(2,2,n-1);
I= eye(2,2);
%  第 n-1层底部反射系数
[Rd(:,:,n)] = ref_single (Vp(n-1),Vs(n-1),rho(n-1),Vp(n),Vs(n),rho(n),seta);

R1(:,:,n-1) = Rd(:,:,n);                      %  （6.38 第1行）   后缀1：下界面    

w = 2*pi*f;
k = w./Vp.*sin(seta);
qeta = sqrt(w.^2.*Vp.^(-2) - k.^2);
yita = sqrt(w.^2.*Vs.^(-2) - k.^2);        %   （6.9）

for j = n-1:-1:2
%          在实际 AVO分析中，角道集数据通常已完成叠前偏移或动校正。因此，为了获得动校正后的角道集合成记录，我们需将相移矩阵(6.32)调整为(6.43)
% E = [exp(-1i*qeta(j).*d(j)) , 0 ; 0 , exp(-1i*yita(j).*d(j))];        %   （6.32）
        E = [exp(-1i*w.*d(j)./Vp(j)) , 0 ; 0 , exp(-1i*w.*d(j)./Vs(j))];      %    (6.43)

R2(:,:,j) = E*R1(:,:,j)*E;             %  第j层顶部反射系数

[Rd(:,:,j),Ru(:,:,j),Td(:,:,j),Tu(:,:,j)] = ref_single (Vp(j-1),Vs(j-1),rho(j-1),Vp(j),Vs(j),rho(j),seta);

R1(:,:,j-1) = Rd(:,:,j) + Tu(:,:,j)*(I - R2(:,:,j)*Ru(:,:,j))^(-1)*R2(:,:,j)*Td(:,:,j);    %  （6.38 第三行）

end

% E = [exp(-1i*qeta(1).*d(1)) , 0 ; 0 , exp(-1i*yita(1).*d(1))];        %   （6.32）
        E = [exp(-1i*w.*d(1)./Vp(1)) , 0 ; 0 , exp(-1i*w.*d(1)./Vs(1))];      %    (6.43)

R2(:,:,1) = E*R1(:,:,1)*E;

[Rd(:,:,1),Ru(:,:,1),Td(:,:,1),Tu(:,:,1)] = ref_single (Vp0,Vs0,rho0,Vp(1),Vs(1),rho(1),seta);

Rm = Rd(:,:,1) + Tu(:,:,1)*(I - R2(:,:,1)*Ru(:,:,1))^(-1)*R2(:,:,1)*Td(:,:,1);    %  （6.38 第三行）


















