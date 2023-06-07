function res = train_fisher_2(X, y)
%输入参数：X为m个因子的n次历史观测值，形状为[n, m] ; y为n次观测的类别（形状为[n, 1]），其值只有0或1两种。
%返回结果：res 为fisher二级判别的结果参数（结构体），包括以下域： 
%res.c  判别系数 (列向量)； 
%res.yc  判别指标 (标量)
%res.meany0  0类样品的y均值(标量)
%res.meany1  1类样品的y均值(标量)
x1 = X(y==1,:);x0 = X(y==0,:);
S1 = (x1-mean(x1,1))'*(x1-mean(x1,1));
S2 = (x0-mean(x0,1))'*(x0-mean(x0,1));
d = mean(x0,1)-mean(x1,1);
c = (S1+S2)\d';
res.c = c/sqrt(sum(c.^2));
res.meany0 = mean(x0,1);
res.meany1 = mean(x1,1);
res.yc = (mean(x0,1)*res.c+mean(x1,1)*res.c)/2;
end

