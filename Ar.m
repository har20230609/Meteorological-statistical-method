function [A,F,p] = Ar(x,n)
% Ar ����Ǳ�׼��������n���Իع�ģ��
% nΪ�Իع�ģ�ͽ���
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


