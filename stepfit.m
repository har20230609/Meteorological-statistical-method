function [A_out,b,bint] = stepfit(A,n,alhpa1,alhpa2)
% stepfit 用于逐步回归
% A矩阵为预报量与因子的离差乘积矩阵/协方差矩阵
% 或A矩阵为标准化的预报量与因子的离差乘积矩阵/协方差矩阵/相关系数矩阵
% b为回归系数
%bint为与检验相关的系数 bint为复相关系数
%若因子数小于等于3 bint(2)为复相关系数的F检验量，bint(3)为1时，回归系数显著区别于0，bint(3)为0时，回归系数不显著
%alhpa1为引入检验参数，alhpa2为剔除检验参数；alhpa1需小于alhpa2；
switch nargin
    case 1
        error('请输入观测次数n')
    case 2
        alhpa1 = 0.05;alhpa2 =0.05;
    case 3
        alhpa2 =0.05;
        if alhpa1 >alhpa2
            disp('未输入剔除检验参数,默认值0.05<%4.2f,引入检验易于剔除检验,引起错误！！',alhpa1)
        end
        if alhpa1>0.3
           error('检查是否输入正确的显著性水平参数：%4.2f>0.3!!',alhpa1)
        end
    case 4
        if alhpa1 >alhpa2
            error('剔除检验参数%4.2f<引入检验参数%4.2f,引入检验易于剔除检验,引起错误！！',alhpa2,alhpa1)
        end
        if alhpa1>0.3||alhpa2>0.3
           error('检查是否输入正确的显著性水平参数：%4.2f>0.3!!',max([alhpa1 alhpa2]))
        end
end
N = length(A);
Syy = A(N,N);
%如果因子数小于等于3，则简单操作
if N<=4
    for i=1:N-1   %对所有因子进行紧凑型变换
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
    
    %逐步回归
else
    
    in0(:) = zeros(N-1,1);                      %in0为该位置是否被引入
    L = length(find(in0==1));               %被引入因子的个数
    for i = 1:3
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %计算方差贡献
        [~,I] = sort(V);                        %对方差贡献进行排序
        in = fliplr(I);                         %in为方差贡献排序
        k = in(find(in0(in)==0,1));             %找到最大方差对应的位置
        in0(k) = 1;                             %引入该位置的因子
        A = cpctf(A,k);                         %对矩阵进行引入变换
        L = length(find(in0==1));               %被引入因子的个数
        F = V(k)/( A(N,N)/(n-L-1) );            %计算方差贡献的检验量
        if F<finv(1-alhpa1,1,n-L-1)
            in0(k) = 0;                         %若引入的因子无法通过引入检验，再剔除该因子
            A = cpctf(A,k);
            L = length(find(in0==1));               %被引入因子的个数
            break
        end
    end
    %最开始引入完成，进行第一次剔除检验
    for j = 1:L-1
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %计算方差贡献
        [~,in] = sort(V);                        %对方差贡献进行排序(倒序)
        k = in(find(in0(in)==1,1));             %找到最小方差对应的位置
        l = length(find(in0==1));               %被引入因子的个数
        F = V(k)/( A(N,N)/(n-l-1) );            %计算方差贡献的检验量
        if F<finv(1-alhpa2,1,n-l-1)
            A = cpctf(A,k);                     %若剔除检验通过
            in0(k) =0;                          %剔除该因子
        else
            break
        end
    end
    
    for i = 4:N-1
        V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %计算方差贡献
        [~,I] = sort(V);                        %对方差贡献进行排序
        in = fliplr(I);                         %in为方差贡献排序
        k = in(find(in0(in)==0,1));             %找到最大方差对应的位置
        in0(k) = 1;                             %引入该位置的因子
        A = cpctf(A,k);                         %对矩阵进行引入变换
        L = length(find(in0==1));               %被引入因子的个数
        F = V(k)/( A(N,N)/(n-L-1) );            %计算方差贡献的检验量
        %进行引入检验
        if F<finv(1-alhpa1,1,n-L-1)
            in0(k) = 0;                         %若引入的因子无法通过引入检验，再剔除该因子
            A = cpctf(A,k);
            break                               %无法再引入因子，退出循环
        else
            %进行剔除检验
            for j = 1:L-1                       %除了新引入的因子其他的因子均需要进行剔除检验
                V(:) = A(1:N-1,N).^2./diag(A(1:N-1,1:N-1));                  %计算方差贡献
                [~,in] = sort(V);                        %对方差贡献进行排序(倒序)
                k = in(find(in0(in)==1,1));             %找到最小方差对应的位置
                l = length(find(in0==1));               %被引入因子的个数
                F = V(k)/( A(N,N)/(n-l-1) );            %计算方差贡献的检验量
                if F<finv(1-alhpa2,1,n-l-1)
                    A = cpctf(A,k);                     %若剔除检验通过
                    in0(k) =0;                          %剔除该因子
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
% 紧凑型变换
function A_out = cpctf(A,n)
% cpctf表示compact transform 用于紧凑型变换
% A为需变换的方阵  n表示变换的位置
switch nargin
    case 1
        if length(A)>1
            error('未输入变换位置');
        else
            error('未输入待变换矩阵');
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


