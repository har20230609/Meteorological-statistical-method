function res = train_bayes(X,y)
% 输入参数：
% X为m个因子的n次历史观测值，形状为[n, m] ;
% y为n次观测的类别（形状为[n, 1]），其值为正整数，取值范围为[1, G], 共有G种类别。
%（为降低难度，约定：n次观测中，从1到G每个类别的资料都存在，建议在程序中对此进行检查，如果不符合该条件，则程序中断报错）。
% p为列向量，形状[G, 1], 表示G个类别的先验概率(任取一样品属于第g类的概率，g=1,2,…,G)。
% 返回结果：res为贝叶斯判别模型的结果参数（结构体类型），包括以下域：
% res.p  形状[G, 1]，各类别的先验概率（如果输入参数p为非空，则与p的值相同）；
% res.c  形状[m, G], 各类别的判别系数c
% res.c0 形状[G,1], 各类别的判别系数c0
V = zeros(size(X,2));
for i=1:max(y)
    p(i)=length(find(y==i))/length(y);
    V = (X(y==i,:)-mean(X(y==i,:),1))'*(X(y==i,:)-mean(X(y==i,:),1))+V;
end
V = V/(size(X,1)-max(y));
for i = 1:max(y)
    c(:,i) = V\(mean(X(y==i,:),1))';
    c0(i) = - mean(X(y==i,:),1)*c(:,i)/2;
end
res.p = p';
res.c = c;
res.c0 = c0';
res.v = V;
end

