function res = train_fisher_2(X, y)
%���������XΪm�����ӵ�n����ʷ�۲�ֵ����״Ϊ[n, m] ; yΪn�ι۲�������״Ϊ[n, 1]������ֵֻ��0��1���֡�
%���ؽ����res Ϊfisher�����б�Ľ���������ṹ�壩������������ 
%res.c  �б�ϵ�� (������)�� 
%res.yc  �б�ָ�� (����)
%res.meany0  0����Ʒ��y��ֵ(����)
%res.meany1  1����Ʒ��y��ֵ(����)
x1 = X(y==1,:);x0 = X(y==0,:);
S1 = (x1-mean(x1,1))'*(x1-mean(x1,1));
S2 = (x0-mean(x0,1))'*(x0-mean(x0,1));
d = mean(x0,1)-mean(x1,1);
c = (S1+S2)\d';
res.c = c/sqrt(sum(c.^2));
res.meany0 = mean(x0,1);
res.meany1 = mean(x1,1);
res.yc = (mean(x0,1)*res.c+mean(x1,1)*res.c)/2;
end

