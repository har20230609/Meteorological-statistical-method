function res = k_means(x,k,sampleDist,initialGroup)
%x��״Ϊn��m�У���ʾn����Ʒ��m��ָ�ꡣ
%k����������ʾ�������
%�������sampleDis,
%�ַ��ͣ���ʾ��Ʒ���룬�ṩ2��ѡ�񣺡�euclidean��(ŷ�Ͼ���)��absolute��(���Ծ���),Ĭ��Ϊŷ�Ͼ��룻
%initialGroup����״Ϊn��1�У�������ʼ���۵㼰�����ȡֵ��ΧΪ0-k��0��ʾ�ǳ�ʼ���۵�(�������)����Ʒ��
%Ҫ��ÿ�������Ҫ�г�ʼ���۵㣬������initialGroup�б�Ȼ����1-k��k����ֵ��
%res: ��״Ϊn��1�У���ʾÿ����Ʒ������������𣬸�Ԫ��ȡֵ��ΧΪ1-k.
n = size(x,1);m = size(x,2);
switch nargin
    case 1
        error('������kֵ')
    case 2
        sampleDist = 'euclidean';
        g1 = randi([0 k],n,1);
    case 3
        g1 = randi([0 k],n,1);
    case 4
        g1 = initialGroup;
end
for i =1:k
    gd(i,1:m) = nanmean(x(g1==i,1:m),1);
end
g0 = zeros(n,1);
while ~isequal(g0,g1)
    g0 = g1;
    for i = 1:n
        g1(i) = groupin(x(i,:),gd,sampleDist);
    end
    for i = 1:k
        gd(i,1:m) = nanmean(x(g1==i,:),1);
    end
end
res = g1;
end

function gin = groupin(d,gd,sampleDist)
switch sampleDist
    case 'euclidean'
        for i = 1:size(gd,1)
            dt(i) = sqrt(sum((d-gd(i,:)).^2));
        end
    case 'absolute'
        for i = 1:size(gd,1)
            dt(i) = sum(abs(d-gd(i,:)));
        end
end
gin = find(abs(dt-min(dt))<1e-3,1);
end
