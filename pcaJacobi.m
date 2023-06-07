function [coeff,score,latent,explained,cum_explained] = pcaJacobi(data,dim)
%pcaJacobi 是用Jacobi法来求不包含NaN数据的矩阵的特征值及特征向量的函数
%  data为m×n的初始矩阵，m为空间格点数，n为时间个数,l为模态数
%输出主成分的系数向量coeff(m*l),主成分score(l*n),特征值latent,方差贡献explained
if any(isnan(data))
    error('数据存在缺省值,无法使用该函数')
end
m = size(data,1);n = size(data,2);
switch nargin
    case 1
        dim = 0;
end
switch dim
    case 0
        data = data-ones(m,n).*mean(data,2);
    case 1
        data = zscore(data,0,2);
end
if m>n
    S = data'*data/(m-1);
    [latent,v] = Jacobi_eigenvalue(S);
%     disp([size(v,1) size(v,2)])
%     disp([size(data,1) size(data,2)])
%     disp([size(sqrt((m-1)*latent(latent~=0)),1) size(sqrt((m-1)*latent(latent~=0)),2)])
    v = data*v(:,latent~=0)./(sqrt(abs((m-1)*latent(latent~=0))))';
    latent = (m-1)/(n-1)*latent;
else
    S = data*data'/(n-1);
    [latent,v] = Jacobi_eigenvalue(S);
end
coeff = v;
explained = latent/sum(latent)*100;
score = coeff'*data;
cum_explained(1) = explained(1);
for i =2:length(explained)
cum_explained(i) =cum_explained(i-1)+explained(i);
end
end

