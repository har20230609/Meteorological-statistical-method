function x_pass = f_pass(x,Tmin,Tmax)
% f_passΪ��ͨ�˲���,�����˹�ָ��Ƶ�ʵ��ź�
% xΪ��ʼ�źš�
%TminΪָ������ĳһ����(���������ݣ�interΪ12��ʾ�˹����ڴ���1����ź�)
%TmaxΪָ��С��ĳһ���ڣ��
% x_passΪ�˹����ź�,TΪ��ͼʱ��Ӧ��ʱ��
n = length(x);
fmin = 1/Tmax;fmax = 1/Tmin;
f = 0:1/(2*n):1/2;
H = ones(1,length(f));
H(f<=fmin|f>=fmax)=0;
h(1) = 1/(2*n)*( H(1)+2*sum( H(2:n)) + H(n+1) );
for i = 1:n
 h(i+1) = 1/(2*n)*( H(1) + 2*sum( H(2:n).*cos(2*pi*i*f(2:n)) ) + H(n+1)*cos(2*pi*i*f(n+1)) );
end
x_pass = zeros(n,1);
for i = 1:n
    c = min(abs(n-i),abs(i-1));
    for j =-c:c
    x_pass(i) = x_pass(i)+h(abs(j)+1)*x(i+j);
    end
end
end

