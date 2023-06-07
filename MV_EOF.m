function [coeff_data,score,latent,explained,cum_explained,lct,s] = MV_EOF(data,lat,lon,latlim,lonlim,dim)
%MV_EOF ����Jacobi�����󲻰���NaN���ݵĶ�������������ֵ�����������ĺ���
%  dataΪf��1�ĳ�ʼ����fΪ�������������е�����ΪMf��n�ľ���
% lat��lonΪnc�����ж�ȡ���ľ�γ�����ݣ�f��m��f��n
% latlim��lonlimΪ����ľ�γ�ȷ�Χ��f��2��f��2
% MfΪ��f�������Ŀռ�������nΪʱ�����
%������ɷֵ�ϵ������coeff������ֵlatent,�����explained
%dim Ĭ��Ϊ0����dimΪ0ʱ�����׼�����coeff_data,��dimΪ1ʱ������Է����coeff_data
if nargin <5
    error('�������㹻����Ϣ')
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

