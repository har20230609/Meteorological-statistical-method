function [gridd,latt,lonn]= getgrid(data,lat,lon,latlim,lonlim,dim)
%getgrid 是获取nc文件中目标格点数据的函数
%   data为三维数据，m×n×t，第一维为纬度,第二维为经度,第三维为时间
% lat与lon为nc数据中读取出的经纬度数据，m×1与n×1
% latlim与lonlim为所需的经纬度范围，2×1与2×1
% dim为0时输出原变量，dim为1时输出距平变量，dim为2时输出标准化变量,dim默认为0
if nargin <5
    error('请输入足够的信息')
elseif nargin==5
    dim = 0;
end
if size(data,1)==length(lon) && size(data,2)==length(lat) && length(lat)~=length(lon)
    data = permute(data,[2,1,3]);
end
gridd = data(lat>=min(latlim)&lat<=max(latlim),lon>=min(lonlim)&lon<=max(lonlim),:);
switch dim
    case 0 
    case 1
        gridd = gridd - nanmean(gridd,3);
    case 2 
        gridd = (gridd-nanmean(gridd,3))./nanstd(gridd,0,3);
end
lonn = lon(lon>=min(lonlim)&lon<=max(lonlim));
latt = lat(lat>=min(latlim)&lat<=max(latlim));
end

