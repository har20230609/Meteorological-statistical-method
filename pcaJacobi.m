function [coeff,score,latent,explained,cum_explained] = pcaJacobi(data,dim)
%pcaJacobi ����Jacobi�����󲻰���NaN���ݵľ��������ֵ�����������ĺ���
%  dataΪm��n�ĳ�ʼ����mΪ�ռ�������nΪʱ�����,lΪģ̬��
%������ɷֵ�ϵ������coeff(m*l),���ɷ�score(l*n),����ֵlatent,�����explained
if any(isnan(data))
    error('���ݴ���ȱʡֵ,�޷�ʹ�øú���')
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

