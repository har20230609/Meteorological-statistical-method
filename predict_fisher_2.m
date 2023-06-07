function class = predict_fisher_2(res, x)
%输入参数：res为fisher二级判别的结果参数(结构体，各域如上所述)；x为m个因子的L次观测（列向量, 形状为[L, m]）
%返回结果：class：标量，表示x的类别，其值为0或1.
%res.c  判别系数 (列向量)；
%res.yc  判别指标 (标量)
%res.meany0  0类样品的y均值(标量)
%res.meany1  1类样品的y均值(标量)
d1 = abs(x*res.c-res.meany1*res.c);
d0 = abs(x*res.c-res.meany0*res.c);
class = zeros(size(x,1),1);
class(d1<d0) = 1;
end

