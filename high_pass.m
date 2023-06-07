function [x_pass,T] = high_pass(x,Tmin,sig)
% low_pass为高通滤波器,用于滤过指定高频率的信号
% x为初始信号，inter为指定大于某一周期(对于月数据，inter为12表示滤过周期大于1年的信号)
% %%本函数的指标为响应函数|H(f)|^2<sig,并且使用二项式系数滑动
% x_pass为滤过的信号,T为画图时对应的时间
n = length(x);
q = find((2*sin(pi/Tmin)).^(1:10^4)>sig,1);
x_pass = zeros(n-q,1);
disp(n-q)
disp(sprintf('使用%d阶差分',q))
for i =q+1:n
    for j = 0:q
        x_pass(i-q) = x_pass(i-q)+(-1)^(j)*nchoosek(q,j)*x(i-j);
    end
end
T = q+1:n;
end