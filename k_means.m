function res = k_means(x,k,sampleDist,initialGroup)
%x形状为n行m列，表示n个样品，m个指标。
%k：整数，表示类别数。
%输入参数sampleDis,
%字符型，表示样品距离，提供2种选择：”euclidean”(欧氏距离)或“absolute”(绝对距离),默认为欧氏距离；
%initialGroup：形状为n行1列，表征初始凝聚点及其类别，取值范围为0-k，0表示非初始凝聚点(暂无类别)的样品。
%要求，每个类别都需要有初始凝聚点，即数组initialGroup中必然包含1-k共k个数值。
%res: 形状为n行1列，表示每个样品最终所属的类别，各元素取值范围为1-k.
n = size(x,1);m = size(x,2);
switch nargin
    case 1
        error('请输入k值')
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
