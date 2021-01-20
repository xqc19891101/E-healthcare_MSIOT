function A = ba(m0,m,N,pp)
if m>m0
    disp('�������m���Ϸ�');  %�������
    return;
end
x=100*rand(1,m0);
y=100*rand(1,m0);
switch  pp  %����ĳ�ʼ��
    case 1
        A=zeros(m0);
    case 2
        A=ones(m0);
        for i=1:m0
            A(i,i)=0;
        end
    case 3
        for i=1:m0
            for j=i+1:m0
                p1=rand(1,1);
                if p1>0.5 
                    A(i,j)=1;A(j,i)=1;
                end
            end
        end
    otherwise
        disp('�������pp���Ϸ�');%�������
        return;          
end 
for k=m0+1:N
    M=size(A,1);
    p=zeros(1,M);
    x0=100*rand(1,1);y0=100*rand(1,1);
    x(k)=x0;y(k)=y0;
    if length(find(A==1))==0
        p(:)=1/M;
    else
         for i=1:M
             p(i)=length(find(A(i,:)==1))/length(find(A==1));
         end
    end
    pp=cumsum(p);          %���ۼƸ���
    for i=1:m              %���ô����еĽڵ������ѡ��m���ڵ����¼���Ľڵ�����
        random_data=rand(1,1);
        aa=find(pp>=random_data);jj=aa(1); 
        A(k,jj)=1;A(jj,k)=1;
    end
end
% plot(x,y,'ro','MarkerEdgeColor','g','MarkerFaceColor','r','markersize',8);hold on;
% for i=1:N
%     for j=i+1:N
%         if A(i,j)~=0
%             plot([x(i),x(j)],[y(i),y(j)],'linewidth',1.2);
%             hold on;
%         end
%     end
% end
% axis equal;hold off;

