function [x_pass,T] = low_pass(x,Tmax,sig)
% low_passΪ��ͨ�˲���,�����˹�ָ����Ƶ�ʵ��ź�
% xΪ��ʼ�źţ�interΪָ������ĳһ����(���������ݣ�interΪ12��ʾ�˹����ڴ���1����ź�)
% %%��������ָ��Ϊ��Ӧ����|H(f)|^2<sig,����ʹ�ö���ʽϵ������
% x_passΪ�˹����ź�,TΪ��ͼʱ��Ӧ��ʱ��
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

