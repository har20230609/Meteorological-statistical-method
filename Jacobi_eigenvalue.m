function [eig,v] = Jacobi_eigenvalue(data)
%Jacobi_eigenvalue ����Jacobi�����󲻰���NaN���ݵľ��������ֵ�����������ĺ���
%
if any(isnan(data))
    error('���ݴ���ȱʡֵ,�޷�ʹ�øú���')
end
n = length(data);
vTemp = zeros(n,n,n); vTemp(vTemp==0)= NaN;
v0 = diag(ones(n,1));
k = 0;
while any(abs(data - diag(diag(data)))>1e-4)
    k = k+1;
    vTemp(:,:,k) = diagonalization2(data);
    data = (squeeze(vTemp(:,:,k)))'*data*(squeeze(vTemp(:,:,k)));
end
for l = 1:k
    v0 = v0*squeeze( vTemp(:,:,l) ) ;
end
[eig,I]= sort(diag(data),'descend');    %��������
v = v0(:,I);
end
function v = diagonalization2(data)
%diagonalization2���ڵõ�����Խǻ�ʱ����������
%�ҵ����ǶԽ���Ԫ�ض�Ӧλ��
data0 = data - diag(diag(data));
[i,j] = find(abs(data0)==max(abs(data0),[],[1 2]),1);
%�������ϵ��
v = diag(ones(length(data),1));
mu = (data(j,j)-data(i,i))/2;deta = data(i,j);
w = mu/abs(mu) * deta/sqrt(mu^2+deta^2);
co = sqrt( (1+sqrt(1-w^2))/2 );
si = w/(2*co);
v(i,i)=co;v(j,j)=co;
v(i,j)=si;v(j,i)=-si;
end