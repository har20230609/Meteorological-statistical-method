function [coeff_data,score,latent,explained,cum_explained,lct,s] = MV_EOF(data,lat,lon,latlim,lonlim,dim)
%MV_EOF 是用Jacobi法来求不包含NaN数据的多变量矩阵的特征值及特征向量的函数
%  data为f×1的初始矩阵，f为变量个数，其中的数据为Mf×n的矩阵
% lat与lon为nc数据中读取出的经纬度数据，f×m与f×n
% latlim与lonlim为所需的经纬度范围，f×2与f×2
% Mf为第f个变量的空间格点数，n为时间个数
%输出主成分的系数向量coeff，特征值latent,方差贡献explained
%dim 默认为0；当dim为0时输出标准化后的coeff_data,当dim为1时输出乘以方差的coeff_data
if nargin <5
    error('请输入足够的信息')
end
f = size(data,2);lct =cell(f,1);
t= size(data{1},3);
Data = [];s = [];
for i =1:f
    [data{i},lct{i}.lat,lct{i}.lon] = getgrid(data{i},lat{i},lon{i},latlim(i,:),lonlim(i,:),0);
    x(i) =size(data{i},1);y(i)=size(data{i},2);
    re_data{i} = reshape(data{i},[x(i)*y(i) t]);
    Mf(i) = x(i)*y(i);
    s = [s;nanstd(re_data{i},0,2)];
    re_data{i} = zscore(re_data{i},0,2);
    Data = [Data;re_data{i}];
end
Data(all(isnan(Data),2),:) = 0;
[coef,score,latent,explained,cum_explained] = pcaJacobi(Data);
coef(coef==0) = NaN;
coeff_data = cell(f,1);
switch dim
    case 0
        coeff_data{1} = reshape(coef(1:Mf(1),:),[x(1) y(1) t]);
        for i = 2:f
            coeff_data{i} = reshape(coef(sum(Mf(1:i-1))+1:sum(Mf(1:i)),:),[x(i) y(i) t]);
        end
    case 1
        coeff_data{1} = reshape(coef(1:Mf(1),:).*s(1:Mf(1)),[x(1) y(1) t]);
        for i = 2:f
            coeff_data{i} = reshape(coef(sum(Mf(1:i-1))+1:sum(Mf(1:i)),:).*s(sum(Mf(1:i-1))+1:sum(Mf(1:i))),[x(i) y(i) t]);
        end
end
end

