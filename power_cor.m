function [s,T,s_alf] = power_cor(x,sig,m)
% �������xΪһάʱ������1*n��sigΪ������ˮƽ,mΪ��󳤶�ѡȡϵ��
% �������s, T ��s_alfΪ������ͬ��һά����
% �ֱ����ͬƵ�ʵĹ�����ֵ�����ڡ����ٽ���ֵ(�����������鷨)��
switch nargin
    case 1
        sig = 0.05;%��δ����������ˮƽ����Ĭ��������ˮƽΪ0.05
        m = 0.5;
    case 2
        m = 0.5;
end
n = length(x);
m1 = floor(m*n);
m2 = ceil(m*n);
if mod(m1,2)==0
    m = m1;
elseif mod(m1,2)==0
    m = m2;
else
    m = m*n+1;
end
r = auto_corr(x,m);
%  r=xcorr(x,x,m,'unbiased');%�����ϵ��
%  r=(r(m+2:end))'; 
fi = 2*pi/m*(1:m-1);
s(1) = 1*(1+2*sum(r(1:m-1))+r(m))/m;
for k = 1:m/2
    s(k+1) = ( 1+2*sum(r(1:m-1).*cos(fi'*k))+r(m)*cos(2*pi*k) )/m;
end

s1(1) = 1/2*(s(1)+s(2));
for k = 1:m/2-1
    s1(k+1) = 1/4*s(k)+1/2*s(k+1)+1/4*s(k+2);
end
s1(m/2+1) = 1/2*(s(m/2)+s(m/2+1));
T = [inf m./(1:m/2)];
s = s1;
s_ba = sum(s1)/(m/2+1);
f = (2*n-m/2)/m;
s_alf = s_ba*chi2inv(1-sig,f)/f*(1-r(1)^2)...
    ./(1-2*r(1)*cos(2*pi*(0:m/2)/m)+r(1)^2);
end
function r = auto_corr(x,t,dim)
% auto_corr ������������غ�������Э��������ĺ���
% x Ϊ 1*n�����У�tΪʱ��
% dimΪ1ʱ�������Ϊ����غ���,dimΪ2ʱ�������Ϊ��Э����� dim Ĭ��Ϊ 1
switch nargin
    case 1
        error('�������㹻����Ϣ')
    case 2
        dim = 1;
end
if dim ==1
    x = zscore(x,0);
end
n = length(x);
if t>n-1
    error('ʱ�ʹ������б�����')
end
t = 1:t;
m = length(t);
r = zeros(m,1);
for i =1:m
    for j =1:n-t(i)
        r(i) = r(i)+(x(j)-nanmean(x))*(x(j+t(i))-nanmean(x));
    end
    r(i) = r(i)/(n-t(i)-1);
end
end

