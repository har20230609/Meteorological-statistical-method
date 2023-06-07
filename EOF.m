function [coeff_data,score,latent,explained,cum_explained,lct] = EOF(data0,lat,lon,latlim,lonlim,dim)
% EOF 是利用getgrid (获取nc文件中目标格点数据的函数)
%和pcaJacobi(用Jacobi法来求不包含NaN数据的矩阵的特征值及特征向量的函数)
%来求解气象数据的主成分与空间系数的函数
% data为三维数据，m×n×t，第一维为纬度,第二维为经度,第三维为时间
% lat与lon为nc数据中读取出的经纬度数据，m×1与n×1
% latlim与lonlim为所需的经纬度范围，2×1与2×1
% 输出主成分的系数向量coeff size(latlim)×size(lonlim)×t
% dim为1时输出距平变量，dim为2时输出标准化变量,dim默认为1
%主成分score,特征值latent,方差贡献explained
if nargin <5
    error('请输入足够的信息')
elseif nargin==5
    dim = 1;
end

[data,lct.latt,lct.lonn] = getgrid(data0,lat,lon,latlim,lonlim,dim);
x =size(data,1);y=size(data,2);t=size(data,3);%得到的数据的三个维度数据
re_data = reshape(data,[x*y t]);
coeff_data = zeros(x*y,t);
n = find(~isnan(re_data(:,1)));
re_data(isnan(re_data(:,1)),:) = [];
[coeff,score,latent,explained,cum_explained] = pcaJacobi(re_data);
coeff_data(n,1:t) = coeff;
coeff_data = reshape(coeff_data,[x y t]);
coeff_data(coeff_data==0) = NaN;
end

