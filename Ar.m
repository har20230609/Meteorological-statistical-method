function [A,F,p] = Ar(x,n)
% Ar 求的是标准化变量的n阶自回归模型
% n为自回归模型阶数
m = length(x);
r= auto_corr(x,m-1);
r = r';
s = diag(ones(n,1));
s(1,:) = [1 r(1:n-1)];
for i =2:n-1
    s(i,:) = [r(i-1:-1:1) 1 r(1:n-i)];
end
s(n,:) = [r(n-1:-1:1) 1];
A = s\(r(1:n))';
F = (sum(A.*(r(1:n))','all')/n)/((1-sum(A.*(r(1:n))','all'))/(m-n-1));
p = fcdf(F,n,m-n-1,'upper');
end



function r = auto_corr(x,t,dim)
% auto_corr 是用于求自相关函数（自协方差函数）的函数
% x 为 1*n的序列，t为时滞
% dim为1时函数结果为自相关函数,dim为2时函数结果为自协方差函数 dim 默认为 1
switch nargin
    case 1
        error('请输入足够的信息')
    case 2
        dim = 1;
end
if dim ==1
x = zscore(x,0);
end
n = length(x); 
if t>n-1
    error('时滞大于序列本身长度')
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


