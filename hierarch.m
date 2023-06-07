function res = hierarch(x,sampleDist,groupDist)
%x��״Ϊn��m�У���ʾn����Ʒ��m��ָ�ꡣ
%�������sampleDis, �ַ��ͣ���ʾ��Ʒ���룬�ṩ2��ѡ�񣺡�euclidean��(ŷ�Ͼ���)��absolute��(���Ծ���)
%�������groupDist, �ַ��ͣ���ʾ�����룬�ṩ3��ѡ�񣺡�maxdist��(������)����mindist��(��С), ��meandist��(ƽ������)
%�������res����״Ϊn��p�С�ÿ�б�ʾһ����Ʒ��res��ʾ����Ʒ��ÿ�����������ϵͳ������ܲ���Ϊp+1��
%res(i,j)��ʾ��i����Ʒ�ڵ�j-1�����������, ���е�1��res(:, 1)ȡֵ��Ϊ[1,2,3,4, ��, n]T����ʾ��0������Ʒ��������𣬼�����Ʒ�Գ�һ�࣬������ͬʱҲ����Ʒ��š�
%�����1���������ǰѵ�4�͵�6����Ʒ�ϲ�Ϊ��n+1�࣬��ô��res��2�е�ֵres(:, 2) = [1, 2, 3, n+1, 5, n+1, ��, n]T�� 
%�����2���������ǰѵ�1, 2����Ʒ�ϲ�Ϊ��n+2�࣬��3,4,6����Ʒ�ϲ�Ϊ��n+3�࣬��ô��res(:,3) = [n+2, n+2, n+3, n+3, 5, n+3, ��, n]T������, 
%��������һ����10��������𶼹�Ϊ1�࣬���������Ϊ20�� ��res����11�У� res(:, 11) ����Ԫ��ֵ��Ϊ20��
switch nargin
    case 1
        sampleDist = 'euclidean';groupDist = 'meandist';
    case 2
        groupDist = 'meandist';
end
n = 1;res(:,n) = (1:size(x,1));
while max(res(:,n))~= min(res(:,n))
    m =max(res(:,n));                          %m���ں�����ŵ�����n��ʾ��n�ξ���
    c = unique(res(:,n));                      %c��ʾ��n�ξ���ʱ���е��������
    n = n+1;  
    dtemp = zeros(length(c));dtemp(dtemp==0)=NaN;
    for i = 1:length(c)
        for j = i+1:length(c)
            dtemp(i,j) = dist(x(res(:,end)==c(i),:),x(res(:,end)==c(j),:),sampleDist,groupDist);
        end
    end
    I = find(dtemp==min(dtemp,[],[1 2]));      %�ҵ���С�����λ��
    a =  ceil(I/size(dtemp,1));b = I - size(dtemp,1)*(a-1);   %�õ���С������ı��
    I0 = judgein(a,b);                                        %�ϲ���
    res(:,n) = res(:,n-1);
    for i = 1:size(I0,1)
        res(any(res(:,n-1)==(c(I0(i,1:find( ~isnan(I0(i,:)),1,'last') )))',2),n)= m+i;
    end
    
end
end

function dt = dist(g1,g2,sampleDist,groupDist)
%dist ���ڼ��������ڸ���֮��ľ��벢�õ�������ĺ���
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
    case 'meandist' %���ƽ������
        dt = mean(d,'all');
end
end

function I = judgein(a,b)
%judgein �����жϽ��ļ��������С���鵱Ϊһ��
%Iÿһ�ж���ʾ��Ϊһ�������ԭ��ţ�I�м��о�˵����Ҫ��Ϊ���顣
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