function [gridd,latt,lonn]= getgrid(data,lat,lon,latlim,lonlim,dim)
%getgrid �ǻ�ȡnc�ļ���Ŀ�������ݵĺ���
%   dataΪ��ά���ݣ�m��n��t����һάΪγ��,�ڶ�άΪ����,����άΪʱ��
% lat��lonΪnc�����ж�ȡ���ľ�γ�����ݣ�m��1��n��1
% latlim��lonlimΪ����ľ�γ�ȷ�Χ��2��1��2��1
% dimΪ0ʱ���ԭ������dimΪ1ʱ�����ƽ������dimΪ2ʱ�����׼������,dimĬ��Ϊ0
if nargin <5
    error('�������㹻����Ϣ')
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

