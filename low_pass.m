function [x_pass,T] = low_pass(x,Tmax,sig)
% low_pass为低通滤波器,用于滤过指定低频率的信号
% x为初始信号，inter为指定大于某一周期(对于月数据，inter为12表示滤过周期大于1年的信号)
% %%本函数的指标为响应函数|H(f)|^2<sig,并且使用二项式系数滑动
% x_pass为滤过的信号,T为画图时对应的时间
n = length(x);
fc = 1/Tmax;
m = find((cos(pi*fc).^(2:2:2*10^4)<sig),2);
Tmax = m(mod(m,2)~=0) -1;
x_pass = zeros(n-Tmax,1);
for j = 0:Tmax
c(j+1) = nchoosek(Tmax,j)/2^Tmax;
end
for i = 1:n-Tmax
    for j = 0:Tmax
    x_pass(i) = x_pass(i)+x(i+j)*c(j+1);
    end
end
T = Tmax/2+1:n-Tmax/2;
disp(m)
end

