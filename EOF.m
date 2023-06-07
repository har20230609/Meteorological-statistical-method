function [coeff_data,score,latent,explained,cum_explained,lct] = EOF(data0,lat,lon,latlim,lonlim,dim)
% EOF ������getgrid (��ȡnc�ļ���Ŀ�������ݵĺ���)
%��pcaJacobi(��Jacobi�����󲻰���NaN���ݵľ��������ֵ�����������ĺ���)
%������������ݵ����ɷ���ռ�ϵ���ĺ���
% dataΪ��ά���ݣ�m��n��t����һάΪγ��,�ڶ�άΪ����,����άΪʱ��
% lat��lonΪnc�����ж�ȡ���ľ�γ�����ݣ�m��1��n��1
% latlim��lonlimΪ����ľ�γ�ȷ�Χ��2��1��2��1
% ������ɷֵ�ϵ������coeff size(latlim)��size(lonlim)��t
% dimΪ1ʱ�����ƽ������dimΪ2ʱ�����׼������,dimĬ��Ϊ1
%���ɷ�score,����ֵlatent,�����explained
if nargin <5
    error('�������㹻����Ϣ')
elseif nargin==5
    dim = 1;
end

[data,lct.latt,lct.lonn] = getgrid(data0,lat,lon,latlim,lonlim,dim);
x =size(data,1);y=size(data,2);t=size(data,3);%�õ������ݵ�����ά������
re_data = reshape(data,[x*y t]);
coeff_data = zeros(x*y,t);
n = find(~isnan(re_data(:,1)));
re_data(isnan(re_data(:,1)),:) = [];
[coeff,score,latent,explained,cum_explained] = pcaJacobi(re_data);
coeff_data(n,1:t) = coeff;
coeff_data = reshape(coeff_data,[x y t]);
coeff_data(coeff_data==0) = NaN;
end

