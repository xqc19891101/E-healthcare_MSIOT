clear;clc
%����BA����
global W
global beita
global mu
A = ba(10,10,500,3);
%��ڵ�Ķ�
userNumber = size(A,1);
nodeDegree = zeros(userNumber,1);
for i = 1:userNumber
    nodeDegree(i) = length(find(A(i,:) == 1));
end
%�ȵķֲ�
uniqueDegree = unique(nodeDegree);
degreeDistribution = zeros(size(uniqueDegree ,1),1);
for d = 1:size(uniqueDegree ,1)
    degreeDistribution(d) = length(find(nodeDegree == uniqueDegree(d))) / userNumber;
end
%%% ����Ӵ�����%%%%
T = 2000; %�Ӵ���������ܵ�ʱ��
flagContact = zeros(userNumber,1); %%�ڵ��ڵ�ǰʱ���Ƿ�Ӵ�
contactEndTime = zeros(userNumber,userNumber); %�ڵ�ÿ�νӴ�����ֹʱ��
contactDuration =  cell(userNumber,userNumber); %�ڵ�ÿ�εĽӴ�����ʱ��
for t = 1:T 
    [row,col] = find(contactEndTime == t); %ÿ��ʱ���ж���û�нڵ�����Ӵ�
    if ~isempty(row) %���ڽ����Ӵ��Ľڵ�
        for i = 1:size(row,1)
            contactEndTime(row(i),col(i)) = 0; %����
            flagContact(row(i)) = 0;  %����
            flagContact(col(i)) = 0;  %����
        end
    end
    for i = 1:userNumber 
        if flagContact(i) == 1  %��ǰʱ���Ƿ�Ӵ�
            continue
        end
        for j = (i+1):userNumber %��ǰʱ��δ�Ӵ�
            if flagContact(j) == 1 %��ǰʱ�̽Ӵ�
                continue
            end
            if A(i,j) ~= 1 %�����û�������neighbor��ϵ
                continue
            end
            if rand <= 0.3  %�����û�����neighbor��ϵ��С��0.3�Ӵ�
                flagContact(i) = 1;
                flagContact(j) = 1;
                contactEndTime(i,j) = t+floor(unifrnd(1,10));
                indoorFlag = 0; %�Ƿ������ڽӴ�
                if rand <= 0.7
                    indoorFlag = 1;
                end
                contactDuration{i,j}(size(contactDuration{i,j},2)+1)= indoorFlag + (contactEndTime(i,j)-t)/10;
                break
            end
        end
    end
end
relationshipType = [1,2,3,4,5];%1:family relationship 2:friendship 3:neighborship 4:colleague relationshio 5:strange relationship
relationshipCell =  cell(userNumber,userNumber); %�����û�֮��Ĺ�ϵ����
for i =  1:userNumber
    for j = (i+1):userNumber
        familyFlag = 0;
        friendFlag = 0;
        neighbourFlag = 0;
        colleagueFlag = 0;
        strangerFlag = 0;
        if rand <=0.4
            strangerFlag = 1;
        end
        if strangerFlag == 0
            if rand <= 0.5
                familyFlag = 1;
            end
            if rand <= 0.5
                friendFlag = 1;
            end
            if rand <= 0.5
                neighbourFlag = 1;
            end
            if rand <= 0.5
                colleagueFlag = 1;
            end
        end
        relationshipCell{i,j}=[familyFlag,friendFlag,neighbourFlag,colleagueFlag,strangerFlag];
    end
end
%%�����Ȩֵ
W = zeros(userNumber,userNumber); %��Ȩֵ����
V = [1.4,1.3,1.2,1.1,1]'; %��ϵ��������
sumV = sum(V)-V(5); %��ϵ���Ⱦ���֮��
yita0 = 0.5; %����Ӵ���Ⱦϵ��
for i =  1:userNumber
    for j = (i+1):userNumber
        contactNumbers = length(contactDuration{i,j});
        partW = 0;
        for n = 1:contactNumbers
            indoorFlag = floor(contactDuration{i,j}(n));
            contactTimelength = (contactDuration{i,j}(n)-indoorFlag)*10;
            if indoorFlag == 1
                 partW = partW + (exp(contactTimelength)-1) * 1;
            else
                 partW = partW + (exp(contactTimelength)-1) * yita0;
            end
        end
        W(i,j) =  (((relationshipCell{i,j}*V)/sumV)*partW);
    end
end
maxW = max(max(W));
for i = 1:userNumber
    for j = (i+1):userNumber
        W(i,j) = W(i,j)/maxW;
        W(j,i) = W(i,j);
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                         ��Ⱦ������ģ��                                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% infectedPercent = 0.1;%��ʼ��ʱ�����д��ڵĸ�Ⱦ����
% infectedNumber = floor(userNumber*infectedPercent);
% randomNumber = randperm(userNumber);
% infectedUsers = randomNumber(1:infectedNumber); %���ѡȡinfectedNumber���ڵ�����Ϊ��ʼ����Ⱦ�ڵ�
% contactUser = zeros(userNumber,1);
% infectedUser = zeros(userNumber,1);
% recoverUser = zeros(userNumber,1);
% healthUser =ones(userNumber,1);
% for i = 1:infectedNumber
%     infectedUser(infectedUsers(i)) = 1;
%     recoverUser(infectedUsers(i)) = 0;
% end   
% beita = 0.08;
% mu = 0.02;
% for t = 1:1000
%     for i = 1:userNumber
%         if healthUser(i) == 1
%             if contactUser(i) == 0
%                 for j = (i+1):userNumber
%                     if contactUser(j) == 0
%                         if infectedUser(j) == 1
%                             if rand < (beita*W(i,j))
%                                 infectedUser(i) 
                            
                        
            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ��Ⱦ����Ⱦģ��                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infectedPercent = 0.1;%��ʼ��ʱ�����д��ڵĸ�Ⱦ����
infectedNumber = floor(userNumber*infectedPercent);
randomNumber = randperm(userNumber);
infectedUsers = randomNumber(1:infectedNumber); %���ѡȡinfectedNumber���ڵ�����Ϊ��ʼ����Ⱦ�ڵ�
initialx =  zeros(1,3*userNumber);
for i = 1:userNumber
    initialx(i) = 1;
    initialx(userNumber+i) = 0;
end
for i = 1:infectedNumber
    initialx(infectedUsers(i)) = 0;
    initialx(userNumber+infectedUsers(i)) = 1;
end
[t,x] = ode45(@infSpread,[0,1000],initialx);
infectedUserNumber = zeros(length(t),1);
for timeStamp = 1:length(t)
    infectedUserNumber(timeStamp) = sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber)));
end
plot(infectedUserNumber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ѡ��ڵ��������                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%1 ��ɫ�㷨�����û���Χ%%%%%%
C = 1:50; %��ʼ������20����ɫ
colorUsed =  zeros(userNumber,1); % ÿ���ڵ�ʹ�õ���ɫ����
saturationDegree = ones(userNumber,1);%�ڵ㱥�Ͷ�
noColorNodes = 1:userNumber; %δȾɫ�ڵ㼯��
while (~isempty(noColorNodes))
    kk = length(noColorNodes)
    maxSaturationNode = find(saturationDegree == max(saturationDegree(noColorNodes)));
    if length(maxSaturationNode) > 1
        nodeDegreeSet = nodeDegree(maxSaturationNode); %�ҵ���󱥺ͶȵĽڵ㼯�ϵĶ�
        colorNode = find(nodeDegree == max(nodeDegreeSet)); %�ҵ����ȵĽڵ�
        colorNode = colorNode(1); %�ҵ���ɫ�ڵ�
    else
        colorNode = maxSaturationNode;  %�ҵ���ɫ�ڵ�
    end
    for c = 1:length(C)
        for j = 1:userNumber
            if A(colorNode,j) > 0 %�����ڵ�����
                if colorUsed(j) == c %%�ڵ�jʹ������ɫc
                    break
                end
            end
        end
        if A(colorNode,j)*colorUsed(j) ~= c %�ڵ�jû���ù���ɫc
            colorUsed(colorNode) = c;
            break
        end
    end
    noColorNodes(find(noColorNodes==colorNode)) = [];
    saturationDegree(colorNode) = 0;
    nodeDegree(colorNode) = -1;
    for i = 1:userNumber %���¼��㱥�Ͷ�
        neighborSet = find(A(i,:) == 1);
        neighborColorSet = colorUsed(neighborSet);
        if saturationDegree(i) ~= 0
            saturationDegree(i) = length(unique(neighborColorSet));
        end
    end  
end
uniColor = unique(colorUsed);
colorNumber = zeros(length(uniColor),1);
for c = 1:length(uniColor)
    colorNumber(c) = length(find(colorUsed == uniColor(c)));
end
maxColor = find(colorNumber == max(colorNumber));
maxColor = uniColor(maxColor(1));
candidateSet = find(colorUsed == maxColor); %%ͨ����ɫ����ѡ���˽ڵ�
%�㷨2����ڵ���Ҫ�ԣ�ѡ�����ױ���Ⱦ�����״�Ⱦ�����˵Ľڵ�
%1)����ڵ�ĸ�Ⱦ����
%�����·�����ߵ����·��������Ȩֵ�ĵ���
W_ = W;
for i = 1:userNumber
    for j = (i+1):userNumber
        if W_(i,j) > 0 
            W_(i,j) = 1/W_(i,j);
            W_(j,i) = W_(i,j);
        end
    end
end
infectedScore = zeros(length(candidateSet),1)%candidate�ı���Ⱦֵ
%��ÿ����ȾԴ���д���
for i = 1:length(infectedUsers) %��ȾԴ
     i
    mark =DJshortpath(W_,infectedUsers(i)); %��ȾԴ��candidate�����·��
    for j = 1:length(candidateSet) %%candidate�еĽڵ�
        j
         [number,Spth] = findPath(mark,infectedUsers(i),candidateSet(j));     
        shortPathNumber = length(Spth); %���·������
        if shortPathNumber == 1 %%���·������Ϊ1
            infected = 1; %%��Ⱦ����
            for k = 1:(length(Spth{1})-1) %���·���ڵ� 
                infected = infected * W(Spth{1}(k),Spth{1}(k+1));
            end
            infectedScore(j)  =  infectedScore(j)+infected; 
        elseif shortPathNumber > 1
            for n = 1:shortPathNumber
                infected = 1;
                for k = 1:(length(Spth{n})-1) %���·���ڵ�
                    infected = infected * W(Spth{n}(k),Spth{n}(k+1));
                end
                 infectedScore(j)  =  infectedScore(j)+infected; 
            end
        end       
    end
end
%2)���Ⱦֵ
infectingScore = zeros(length(candidateSet),1);
for i = 1:length(candidateSet)
    shortPathHopSet = cell(userNumber-1,1);
     mark =DJshortpath(W_,candidateSet(i)); %��ȾԴ��candidate�����·��
     for j = 1:userNumber
         dis=mark{j,2};
         if dis == inf
             break
         end
         [number,Spth] = findPath(mark,candidateSet(i),j);
          shortPathNumber = length(Spth); %���·������
          if shortPathNumber == 1 %%���·������Ϊ1
              hopNumber = length(Spth{1})-1;
              shortPathHopSet{hopNumber}{length(shortPathHopSet{hopNumber})+1} = Spth{1};
          elseif shortPathNumber > 1
               for n = 1:shortPathNumber
                   hopNumber = length(Spth{n})-1;
                   shortPathHopSet{hopNumber}{length(shortPathHopSet{hopNumber})+1} = Spth{n};
               end
          end
     end  
     for h = 1:length(shortPathHopSet)
         pathNumber = length(shortPathHopSet{h});
         if pathNumber > 0
             for n = 1:pathNumber
                 infecting = 1;
                 for k = 1:h
                     infecting = infecting*W(shortPathHopSet{h}{n}(k),shortPathHopSet{h}{n}(k+1));
                 end
                 infectingScore(i) = infectingScore(i) + infecting;
             end
         end
     end               
end
combinedScore = zeros(length(candidateSet),1);
for i = 1:length(candidateSet)
    combinedScore(i) = (infectedScore(i)/max(infectedScore))*(infectingScore(i)/max(infectingScore));
    if ~isempty(find(infectedUsers == candidateSet(i)))
        combinedScore(i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beita = 0.020;
mu = 0.02;
k = 1;
lastInfectedUserNumber = [];
for immuneNumber = 0:10:50
    combinedScoreTempt = combinedScore;
    sortCombinedScoreTempt = sort(combinedScoreTempt,'descend');
    immuneNodeSet = zeros(immuneNumber,1);
    for i = 1:immuneNumber
        immuneNodePlace = find(combinedScore == sortCombinedScoreTempt(i));
        for j = 1:length(immuneNodePlace)
            immuneNodeTempt = candidateSet(immuneNodePlace(j));
            if isempty(find(immuneNodeSet == immuneNodeTempt))
                immuneNodeSet(i) = immuneNodeTempt;
                break
            end
        end
    end
    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)))
    k=k+1;
    
end
plot(lastInfectedUserNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%random scheme %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
lastInfectedUserNumber = [];
for immuneNumber = 0:7:35
    immuneNodeSet = zeros(immuneNumber,1);
    randomNumber = randperm(userNumber);
    for i = 1:length(infectedUsers)
        randomNumber(find(randomNumber == infectedUsers(i))) = [];
    end
     immuneNodeSet = randomNumber(1:immuneNumber); %���ѡȡinfectedNumber���ڵ�����Ϊ��ʼ����Ⱦ�ڵ� 

    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)),'r')
    k = k+1;
    
end
plot(lastInfectedUserNumber,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%degree scheme %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
lastInfectedUserNumber = [];
for immuneNumber =  0:7:35
    immuneNodeSet = zeros(immuneNumber,1);
    userNumber = size(A,1);
    nodeDegree = zeros(userNumber,1);
    for i = 1:userNumber
        nodeDegree(i) = length(find(A(i,:) == 1));
    end
    combinedScoreTempt = nodeDegree;
    combinedScoreTempt(infectedUsers) = 0;
    sortCombinedScoreTempt = sort(combinedScoreTempt,'descend');
    immuneNodeSet = zeros(immuneNumber,1);
    for i = 1:immuneNumber
        immuneNodePlace = find(nodeDegree == sortCombinedScoreTempt(i));
        for j = 1:length(immuneNodePlace)
            immuneNodeTempt = immuneNodePlace(j);
            if isempty(find(immuneNodeSet == immuneNodeTempt))
                immuneNodeSet(i) = immuneNodeTempt;
                break
            end
        end
    end
    

    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)),'r')
    k = k+1;
    
end
plot(lastInfectedUserNumber,'g')
              
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Ƚ�beita��Ӱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%beita = 0.020;
mu = 0.02;
k = 1;
lastInfectedUserNumber = [];
for beita = 0.010:0.002:0.02
    immuneNumber = 30
    combinedScoreTempt = combinedScore;
    sortCombinedScoreTempt = sort(combinedScoreTempt,'descend');
    immuneNodeSet = zeros(immuneNumber,1);
    for i = 1:immuneNumber
        immuneNodePlace = find(combinedScore == sortCombinedScoreTempt(i));
        for j = 1:length(immuneNodePlace)
            immuneNodeTempt = candidateSet(immuneNodePlace(j));
            if isempty(find(immuneNodeSet == immuneNodeTempt))
                immuneNodeSet(i) = immuneNodeTempt;
                break
            end
        end
    end
    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)))
    k=k+1;
    
end
plot(lastInfectedUserNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%random scheme %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
lastInfectedUserNumber = [];
for beita = 0.010:0.002:0.02
    immuneNumber = 20
    immuneNodeSet = zeros(immuneNumber,1);
    randomNumber = randperm(userNumber);
    for i = 1:length(infectedUsers)
        randomNumber(find(randomNumber == infectedUsers(i))) = [];
    end
     immuneNodeSet = randomNumber(1:immuneNumber); %���ѡȡinfectedNumber���ڵ�����Ϊ��ʼ����Ⱦ�ڵ� 

    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)),'r')
    k = k+1;
    
end
plot(lastInfectedUserNumber,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%degree scheme %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
lastInfectedUserNumber = [];
for beita = 0.010:0.002:0.02
    immuneNumber =20
    immuneNodeSet = zeros(immuneNumber,1);
    userNumber = size(A,1);
    nodeDegree = zeros(userNumber,1);
    for i = 1:userNumber
        nodeDegree(i) = length(find(A(i,:) == 1));
    end
    combinedScoreTempt = nodeDegree;
    combinedScoreTempt(infectedUsers) = 0;
    sortCombinedScoreTempt = sort(combinedScoreTempt,'descend');
    immuneNodeSet = zeros(immuneNumber,1);
    for i = 1:immuneNumber
        immuneNodePlace = find(nodeDegree == sortCombinedScoreTempt(i));
        for j = 1:length(immuneNodePlace)
            immuneNodeTempt = immuneNodePlace(j);
            if isempty(find(immuneNodeSet == immuneNodeTempt))
                immuneNodeSet(i) = immuneNodeTempt;
                break
            end
        end
    end
    

    initialx =  zeros(1,3*userNumber);
    for i = 1:userNumber
        initialx(i) = 1;
        initialx(userNumber+i) = 0;
    end
    for i = 1:infectedNumber
        initialx(infectedUsers(i)) = 0;
        initialx(userNumber+infectedUsers(i)) = 1;
    end
    for i = 1:immuneNumber
        initialx(immuneNodeSet(i)) = 0;
        initialx(userNumber+immuneNodeSet(i)) = 0;
    end
    [t,x] = ode45(@infSpread,[0,1000],initialx);
    infectedUserNumber = zeros(length(t),1);
    for timeStamp = 1:length(t)
        infectedUserNumber(timeStamp) = (sum(x(timeStamp,(userNumber+1):(2*userNumber)))+sum(x(timeStamp,(2*userNumber+1):(3*userNumber))))/userNumber;
    end
    hold on
    lastInfectedUserNumber(k) = infectedUserNumber(timeStamp) ;
    %plot(t(1:10:length(t)),infectedUserNumber(1:10:length(t)),'r')
    k = k+1;
    
end
plot(lastInfectedUserNumber,'g')
              





        
        

            

    
