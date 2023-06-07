function res = hierarch(x,sampleDist,groupDist)
%x形状为n行m列，表示n个样品，m个指标。
%输入参数sampleDis, 字符型，表示样品距离，提供2种选择：”euclidean”(欧氏距离)或“absolute”(绝对距离)
%输入参数groupDist, 字符型，表示类间距离，提供3种选择：“maxdist”(最大距离)，”mindist”(最小), 或“meandist”(平均距离)
%结果变量res：形状为n行p列。每行表示一个样品，res表示各样品在每步所属的类别。系统聚类的总步数为p+1。
%res(i,j)表示第i个样品在第j-1步所属的类别, 其中第1列res(:, 1)取值恒为[1,2,3,4, …, n]T，表示第0步各样品所属的类别，即各样品自成一类，类别号码同时也是样品序号。
%如果第1步聚类结果是把第4和第6个样品合并为第n+1类，那么，res第2列的值res(:, 2) = [1, 2, 3, n+1, 5, n+1, …, n]T； 
%如果第2步聚类结果是把第1, 2个样品合并为第n+2类，第3,4,6个样品合并为第n+3类，那么，res(:,3) = [n+2, n+2, n+3, n+3, 5, n+3, …, n]T，……, 
%如果在最后一步第10步所有类别都归为1类，新类别命名为20， 则res共有11列， res(:, 11) 各行元素值均为20。
switch nargin
    case 1
        sampleDist = 'euclidean';groupDist = 'meandist';
    case 2
        groupDist = 'meandist';
end
n = 1;res(:,n) = (1:size(x,1));
while max(res(:,n))~= min(res(:,n))
    m =max(res(:,n));                          %m用于后续组号调整，n表示第n次聚类
    c = unique(res(:,n));                      %c表示第n次聚类时所有的组号向量
    n = n+1;  
    dtemp = zeros(length(c));dtemp(dtemp==0)=NaN;
    for i = 1:length(c)
        for j = i+1:length(c)
            dtemp(i,j) = dist(x(res(:,end)==c(i),:),x(res(:,end)==c(j),:),sampleDist,groupDist);
        end
    end
    I = find(dtemp==min(dtemp,[],[1 2]));      %找到最小距离的位置
    a =  ceil(I/size(dtemp,1));b = I - size(dtemp,1)*(a-1);   %得到最小距离组的编号
    I0 = judgein(a,b);                                        %合并组
    res(:,n) = res(:,n-1);
    for i = 1:size(I0,1)
        res(any(res(:,n-1)==(c(I0(i,1:find( ~isnan(I0(i,:)),1,'last') )))',2),n)= m+i;
    end
    
end
end

function dt = dist(g1,g2,sampleDist,groupDist)
%dist 用于计算两组内各点之间的距离并得到组间距离的函数
d = zeros(length(g1),length(g2));
for i = 1:size(g1,1)
    for j = 1:size(g2,1)
        switch sampleDist
            case 'euclidean'
                d(i,j) = sqrt(sum((g1(i,:)-g2(j,:)).^2));
            case 'absolute'
                d(i,j) = sum(abs(g1(i,:)-g2(j,:)));
        end
    end
end
switch groupDist
    case 'maxdist'
        dt = max(d,[],[1 2]);
    case 'mindist'
        dt = min(d,[],[1 2]);
    case 'meandist' %组间平均距离
        dt = mean(d,'all');
end
end

function I = judgein(a,b)
%judgein 用于判断将哪几类距离最小的组当为一组
%I每一行都表示分为一个新组的原组号，I有几行就说明需要分为几组。
I = zeros(length(a),2*length(a));I(I==0)=NaN;
I(1,1:2) = [a(1) b(1)];
for i = 2:length(a)
    if any(a(i) ==I)
        n = find(any(I==a(i),2));
        m = find(isnan(I(n,:)));
        I(n,m:m+1) = [a(i),b(i)];
    elseif any(b(i) ==I)
        n = find(any(I==b(i),2));
        m = find(isnan(I(n,:)));
        I(n,m:m+1) = [a(i),b(i)];
    else
        I(find(isnan(I),1),1:2) = [a(i),b(i)];
    end
end
I(all(isnan(I),2),:)=[];
end