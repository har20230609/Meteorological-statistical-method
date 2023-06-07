function  [s,T,s_alf] = power_discrt(x,sig)
% 输入参数x为一维时间序列1*n，sig为显著性水平
% 输出参数s, T 和s_alf为长度相同的一维数组
% 分别代表不同频率的功率谱值、周期、和临界谱值(按红噪声检验法)。
switch nargin
    case 1
        sig = 0.05;%若未输入显著性水平，则默认显著性水平为0.05
end
n = length(x);
fi = 2*pi/n*(1:n);
for k =1:floor(n/2)
a(k) = 2/n*sum(x.*cos(fi*k));
b(k) = 2/n*sum(x.*sin(fi*k));
end
if mod(n,2)==0  %如果n为偶数，进行调整
a(n/2) = 1/n*sum(x.*cos(fi*n/2));
b(n/2) = 0;
end
s = (a.^2+b.^2)/2;
T = n./(1:floor(n/2));
Fc = finv(1-sig,2,n-3);
s_alf = (2*Fc*var(x)/(2*Fc+n-3))*ones(1,floor(n/2));
end

