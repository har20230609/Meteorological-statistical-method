function [x_pass,T] = high_pass(x,Tmin,sig)
% low_passΪ��ͨ�˲���,�����˹�ָ����Ƶ�ʵ��ź�
% xΪ��ʼ�źţ�interΪָ������ĳһ����(���������ݣ�interΪ12��ʾ�˹����ڴ���1����ź�)
% %%��������ָ��Ϊ��Ӧ����|H(f)|^2<sig,����ʹ�ö���ʽϵ������
% x_passΪ�˹����ź�,TΪ��ͼʱ��Ӧ��ʱ��
n = length(x);
q = find((2*sin(pi/Tmin)).^(1:10^4)>sig,1);
x_pass = zeros(n-q,1);
disp(n-q)
disp(sprintf('ʹ��%d�ײ��',q))
for i =q+1:n
    for j = 0:q
        x_pass(i-q) = x_pass(i-q)+(-1)^(j)*nchoosek(q,j)*x(i-j);
    end
end
T = q+1:n;
end