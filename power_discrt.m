function  [s,T,s_alf] = power_discrt(x,sig)
% �������xΪһάʱ������1*n��sigΪ������ˮƽ
% �������s, T ��s_alfΪ������ͬ��һά����
% �ֱ����ͬƵ�ʵĹ�����ֵ�����ڡ����ٽ���ֵ(�����������鷨)��
switch nargin
    case 1
        sig = 0.05;%��δ����������ˮƽ����Ĭ��������ˮƽΪ0.05
end
n = length(x);
fi = 2*pi/n*(1:n);
for k =1:floor(n/2)
a(k) = 2/n*sum(x.*cos(fi*k));
b(k) = 2/n*sum(x.*sin(fi*k));
end
if mod(n,2)==0  %���nΪż�������е���
a(n/2) = 1/n*sum(x.*cos(fi*n/2));
b(n/2) = 0;
end
s = (a.^2+b.^2)/2;
T = n./(1:floor(n/2));
Fc = finv(1-sig,2,n-3);
s_alf = (2*Fc*var(x)/(2*Fc+n-3))*ones(1,floor(n/2));
end

