function [A_out,b,bint] = stepfit(A,n,alhpa1,alhpa2)
% stepfit �����𲽻ع�
% A����ΪԤ���������ӵ����˻�����/Э�������
% ��A����Ϊ��׼����Ԥ���������ӵ����˻�����/Э�������/���ϵ������
% bΪ�ع�ϵ��
%bintΪ�������ص�ϵ�� bintΪ�����ϵ��
%��������С�ڵ���3 bint(2)Ϊ�����ϵ����F��������bint(3)Ϊ1ʱ���ع�ϵ������������0��bint(3)Ϊ0ʱ���ع�ϵ��������
%alhpa1Ϊ������������alhpa2Ϊ�޳����������alhpa1��С��alhpa2��
switch nargin
    case 1
        error('������۲����n')
    case 2
        alhpa1 = 0.05;alhpa2 =0.05;
    case 3
        alhpa2 =0.05;
        if alhpa1 >alhpa2
            disp('δ�����޳��������,Ĭ��ֵ0.05<%4.2f,������������޳�����,������󣡣�',alhpa1)
        end
        if alhpa1>0.3
           error('����Ƿ�������ȷ��������ˮƽ������%4.2f>0.3!!',alhpa1)
        end
    case 4
        if alhpa1 >alhpa2
            error('�޳��������%4.2f<����������%4.2f,������������޳�����,������󣡣�',alhpa2,alhpa1)
        end
        if alhpa1>0.3||alhpa2>0.3
           error('����Ƿ�������ȷ��������ˮƽ������%4.2f>0.3!!',max([alhpa1 alhpa2]))
        end
end
N = length(A);
Syy = A(N,N);
%���������С�ڵ���3����򵥲���
if N<=4
    for i=1:N-1   %���������ӽ��н����ͱ任
        A = cpctf(A,i);
    end
    A_out = A;
    bint(1) = 1 - A(N,N)/Syy;
    bint(2) = (bint(1)/(N-1)) /( (1-bint(1))/(n-N) );
    bint(4) = fcdf(bint(2),N-1,n-N,'upper');
    if bint(2)>finv(1-alhpa1,N-1,n-N)
        bint(3) = 1;
    else
        bint(3) = 0;
    end
    b = A_out(1:N-1,N);
    
    %�𲽻ع�
else
    
    in0(:) = zeros(N-1,1);                      %in0Ϊ��λ���Ƿ�����
    L = length(find(in0==1));               %���������ӵĸ���
    for i = 1:3
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %���㷽���
        [~,I] = sort(V);                        %�Է���׽�������
        in = fliplr(I);                         %inΪ���������
        k = in(find(in0(in)==0,1));             %�ҵ���󷽲��Ӧ��λ��
        in0(k) = 1;                             %�����λ�õ�����
        A = cpctf(A,k);                         %�Ծ����������任
        L = length(find(in0==1));               %���������ӵĸ���
        F = V(k)/( A(N,N)/(n-L-1) );            %���㷽��׵ļ�����
        if F<finv(1-alhpa1,1,n-L-1)
            in0(k) = 0;                         %������������޷�ͨ��������飬���޳�������
            A = cpctf(A,k);
            L = length(find(in0==1));               %���������ӵĸ���
            break
        end
    end
    %�ʼ������ɣ����е�һ���޳�����
    for j = 1:L-1
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %���㷽���
        [~,in] = sort(V);                        %�Է���׽�������(����)
        k = in(find(in0(in)==1,1));             %�ҵ���С�����Ӧ��λ��
        l = length(find(in0==1));               %���������ӵĸ���
        F = V(k)/( A(N,N)/(n-l-1) );            %���㷽��׵ļ�����
        if F<finv(1-alhpa2,1,n-l-1)
            A = cpctf(A,k);                     %���޳�����ͨ��
            in0(k) =0;                          %�޳�������
        else
            break
        end
    end
    
    for i = 4:N-1
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %���㷽���
        [~,I] = sort(V);                        %�Է���׽�������
        in = fliplr(I);                         %inΪ���������
        k = in(find(in0(in)==0,1));             %�ҵ���󷽲��Ӧ��λ��
        in0(k) = 1;                             %�����λ�õ�����
        A = cpctf(A,k);                         %�Ծ����������任
        L = length(find(in0==1));               %���������ӵĸ���
        F = V(k)/( A(N,N)/(n-L-1) );            %���㷽��׵ļ�����
        %�����������
        if F<finv(1-alhpa1,1,n-L-1)
            in0(k) = 0;                         %������������޷�ͨ��������飬���޳�������
            A = cpctf(A,k);
            break                               %�޷����������ӣ��˳�ѭ��
        else
            %�����޳�����
            for j = 1:L-1                       %������������������������Ӿ���Ҫ�����޳�����
                V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %���㷽���
                [~,in] = sort(V);                        %�Է���׽�������(����)
                k = in(find(in0(in)==1,1));             %�ҵ���С�����Ӧ��λ��
                l = length(find(in0==1));               %���������ӵĸ���
                F = V(k)/( A(N,N)/(n-l-1) );            %���㷽��׵ļ�����
                if F<finv(1-alhpa2,1,n-l-1)
                    A = cpctf(A,k);                     %���޳�����ͨ��
                    in0(k) =0;                          %�޳�������
                else
                    break
                end
            end
        end
    end
    
    A_out = A;
    b = A_out(1:N-1,N);
    b(in0==0) = 0;
    bint = 1 - A_out(N,N)/Syy;
end
end
% �����ͱ任
function A_out = cpctf(A,n)
% cpctf��ʾcompact transform ���ڽ����ͱ任
% AΪ��任�ķ���  n��ʾ�任��λ��
switch nargin
    case 1
        if length(A)>1
            error('δ����任λ��');
        else
            error('δ������任����');
        end
    case 2
        N = length(A);
        A_out = ones(N,N);
        A_out(n,n) = 1/A(n,n);
        A_out(n,(1:N)~=n) = A(n,(1:N)~=n)/A(n,n);
        A_out((1:N)~=n,n) = -A((1:N)~=n,n)/A(n,n);
        for i =1:N
            for j =1:N
                if i == n || j==n
                    continue
                else
                    A_out(i,j) = A(i,j) - A(n,j)*A(i,n)/A(n,n);
                end
            end
        end
end
end


