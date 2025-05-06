function wigb (a,scal,x,z,amx)
%WIGB: Plot seismic data using wiggles ʹ�ðڶ����Ƶ�������
%���볹��Ū�����������
%  WIGB(a,scal,x,z,amx) 
%
%  IN    a: seismic data
%        scale: multiple data by scale+
%        x: x-axis;
%        z: vertical axis (time or depth)
%	 x and z are vectors with offset and time.
%
%	 If only 'a' is enter, 'scal,x,z,amn,amx' are decided automatically; 
%	 otherwise, 'scal' is a scalar; 'x, z' are vectors for annotation in 
%	 offset and time, amx are the amplitude range.
%
% Author:
% 	Xingong Li, Dec. 1995
% Changes: 
%	Jun11,1997: add amx
% 	May16,1997: updated for v5 - add 'zeros line' to background color
% 	May17,1996: if scal ==0, plot without scaling
% 	Aug6, 1996: if max(tr)==0, plot a line 
%   07-7-13,�붮

if nargin == 0       %Ĭ�������10��,10������ 
    nx=10;nz=10; 
    a = rand(nz,nx)-0.5; %�������-0.5~0.5֮��
end;

[nz,nx]=size(a); %nz��������,nx����

trmx= max(abs(a));   %trmx��¼ÿһ���ľ���ֵ�����ֵ
if (nargin <= 4); amx=mean(trmx);  end;%Ĭ�������,amx��Ϊ�����������ֵ���ֵ��ƽ��ֵ
if (nargin <= 2); x=[1:nx]; z=[1:nz]; end;%Ĭ������µ�x,zֵ
if (nargin <= 1); scal =1; end;%Ĭ������µ����ű���Ϊ1

if nx <= 1; disp(' ERR:PlotWig: nx has to be more than 1');return;end;%����С�ڵ���1��,��ʾ������Ϣ

 % take the average as dx
	dx1 = abs(x(2:nx)-x(1:nx-1));
 	dx = median(dx1);%�ڱ�עx������ʱ,���ݱ�ע��С,������������ݴ���һ��

 dz=z(2)-z(1);
 xmx=max(a(:)); xmn=min(a(:)); %�������ݵ�����ֵ����Сֵ

 if scal == 0; scal=1; end;%���ű��趨Ϊ1
 a = a * dx /amx; %��û��ָ��x,z�������,dx=1
 a = a * scal;

 fprintf(' PlotWig: data range [%f, %f], plotted max %f \n',xmn,xmx,amx);
 %-------------------------------------------------------------------------

% set display range 
x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
z1=min(z)-dz; z2=max(z)+dz;
 
set(gca,'NextPlot','add','Box','on', ...
   'XAxisLocation','bottom', ...
  'XLim', [x1 x2], ...
  'YDir','reverse', ...
  'YLim',[z1 z2]);%�����������һЩ����,���������������ϵ��������Сֵ,x������,��math����.ѧ���Լ���set(gca)
 

	fillcolor = [0 0 0];
	linecolor = [0 0 0];
	linewidth = 0.1;

	z=z'; 	% input as row vector
	zstart=z(1);
	zend  =z(nz);

for i=1:nx,
   
  if trmx(i) ~= 0;    % skip the zero traces
	tr=a(:,i); 	% --- one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% ���Ų�ͬ,˵��ʲô,�й���㰡%�޸���s(1:nz-1) ~= (-s(2:nz)
	npos = length(i1);%������~


	%12/7/97 
	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
	aadd = zeros(size(zadd));
    if size(zadd==npos) 
        disp('zadd==nops');
    end

	zpos = find(tr >0);    %ɾȥ[zpos,vpos]�е�vpos
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
		if tr(1)>0, 	a0=0; z0=1.00;
        else 		a0=0; z0=zadd(1);
		end;
		if tr(nz)>0, 	a1=0; z1=nz; 
        else 		a1=0; z1=max(zadd);
		end;
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		

	zzz = zstart + zz*dz -dz;

	patch( aa+x(i) , zzz,  fillcolor);

	line( 'Color',[1 1 1],  ...
         'Xdata', x(i)+[0 0], 'Ydata',[zstart zend]); % remove zero line

%'LineWidth',linewidth, ...
%12/7/97 'Xdata', x(i)+[0 0], 'Ydata',[z0 z1]*dz);	% remove zero line

	line( 'Color',linecolor,  ...
	 'LineWidth',linewidth, ...
	 'Xdata', tr+x(i), 'Ydata',z);	% negatives line

   else % zeros trace
	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
         'Xdata', [x(i) x(i)], 'Ydata',[zstart zend]);
   end;
end;
% title 'ʱ������'
%colorbar







