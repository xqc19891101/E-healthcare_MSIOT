function mark=DJshortpath(A,st)
% Spth=st(2);
n=length(A);
mark=cell(n,2);
for i = 1:n
    mark{i,1} =  st(1);
    mark{i,2} = inf;
    mark{st(1),2} = 0;
end
P=st(1);%P?????邦??????P㊣那??????.
TP=st(1);%TP?????邦??????T㊣那??㊣???P㊣那??????.
m=0;
while length(P)<n %?迄?????迄??P㊣那??,㊣那???芍??.
    for k=1:length(TP)
        for i=1:n
            if A(TP(k),i)~=0
                if mark{TP(k),2}+A(TP(k),i)<mark{i,2}
                    %????T㊣那??.
                    mark{i,2}=mark{TP(k),2}+A(TP(k),i);
                    mark{i,1}=TP(k);
                elseif mark{TP(k),2}+A(TP(k),i)==mark{i,2}
                    mark{i,1} = [mark{i,1},TP(k)];
                end
            end
        end
    end 
    M=[mark{:,2}];
    m=min(M([mark{:,2}]>m));%?車?迄??T㊣那????℅?????.
    TP=find([mark{:,2}]==m);%?谷??T㊣那????℅???????.
    P=[P TP];%??T㊣那????℅???????T㊣那??????P㊣那??.
end
%?車℅??????∟?∼????.
% d=mark{st(2),2};
% [number,Spth] = findPath(mark,st(1),st(2));

