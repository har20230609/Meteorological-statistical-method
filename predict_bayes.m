function [class,p_post] = predict_bayes (res, x) 
%���������resΪ��Ҷ˹�б�Ľ��������xΪm�����ӵ�L�ι۲⣨������, ��״Ϊ[m,L]��
%���ؽ����
%class����������ʾx����� 
%p_post: ��������״[G,L]��������ʣ���ʾ��Ʒx���ڸ����ĸ��ʣ�Ҳ���������ʣ���
yg = log(res.p)+(res.c)'*x+res.c0;
for i =1:length(res.c0)
p_post(i,:) = 1./sum(exp(yg-yg(i,:)),1);
end
[~,class] = max(p_post,[],1);
end

