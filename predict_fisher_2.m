function class = predict_fisher_2(res, x)
%���������resΪfisher�����б�Ľ������(�ṹ�壬������������)��xΪm�����ӵ�L�ι۲⣨������, ��״Ϊ[L, m]��
%���ؽ����class����������ʾx�������ֵΪ0��1.
%res.c  �б�ϵ�� (������)��
%res.yc  �б�ָ�� (����)
%res.meany0  0����Ʒ��y��ֵ(����)
%res.meany1  1����Ʒ��y��ֵ(����)
d1 = abs(x*res.c-res.meany1*res.c);
d0 = abs(x*res.c-res.meany0*res.c);
class = zeros(size(x,1),1);
class(d1<d0) = 1;
end

