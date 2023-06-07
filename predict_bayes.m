function [class,p_post] = predict_bayes (res, x) 
%输入参数：res为贝叶斯判别的结果参数，x为m个因子的L次观测（列向量, 形状为[m,L]）
%返回结果：
%class：标量，表示x的类别； 
%p_post: 列向量形状[G,L]，后验概率，表示样品x属于各类别的概率（也即条件概率）。
yg = log(res.p)+(res.c)'*x+res.c0;
for i =1:length(res.c0)
p_post(i,:) = 1./sum(exp(yg-yg(i,:)),1);
end
[~,class] = max(p_post,[],1);
end

