clear;clc
%建立BA网络
global W
global beita
global mu
A = ba(10,10,500,3);
%求节点的度
userNumber = size(A,1);
nodeDegree = zeros(userNumber,1);
for i = 1:userNumber
    nodeDegree(i) = length(find(A(i,:) == 1));
end
%度的分布
uniqueDegree = unique(nodeDegree);
degreeDistribution = zeros(size(uniqueDegree ,1),1);
for d = 1:size(uniqueDegree ,1)
    degreeDistribution(d) = length(find(nodeDegree == uniqueDegree(d))) / userNumber;
end
%%% 构造接触数据%%%%
T = 2000; %接触输出构造总的时间
flagContact = zeros(userNumber,1); %%节点在当前时刻是否接触
contactEndTime = zeros(userNumber,userNumber); %节点每次接触的终止时刻
contactDuration =  cell(userNumber,userNumber); %节点每次的接触持续时间
for t = 1:T 
    [row,col] = find(contactEndTime == t); %每个时刻判断有没有节点结束接触
    if ~isempty(row) %存在结束接触的节点
        for i = 1:size(row,1)
            contactEndTime(row(i),col(i)) = 0; %清零
            flagContact(row(i)) = 0;  %清零
            flagContact(col(i)) = 0;  %清零
        end
    end
    for i = 1:userNumber 
        if flagContact(i) == 1  %当前时刻是否接触
            continue
        end
        for j = (i+1):userNumber %当前时刻未接触
            if flagContact(j) == 1 %当前时刻接触
                continue
            end
            if A(i,j) ~= 1 %两个用户不存在neighbor关系
                continue
            end
            if rand <= 0.3  %两个用户存在neighbor关系，小于0.3接触
                flagContact(i) = 1;
                flagContact(j) = 1;
                contactEndTime(i,j) = t+floor(unifrnd(1,10));
                indoorFlag = 0; %是否在室内接触
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
relationshipCell =  cell(userNumber,userNumber); %所有用户之间的关系矩阵
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
%%计算边权值
W = zeros(userNumber,userNumber); %边权值矩阵
V = [1.4,1.3,1.2,1.1,1]'; %关系力度向量
sumV = sum(V)-V(5); %关系力度矩阵之和
yita0 = 0.5; %室外接触感染系数
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
% %                         传染病传播模拟                                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% infectedPercent = 0.1;%初始化时网络中存在的感染比例
% infectedNumber = floor(userNumber*infectedPercent);
% randomNumber = randperm(userNumber);
% infectedUsers = randomNumber(1:infectedNumber); %随机选取infectedNumber个节点来作为初始化感染节点
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
%                         传染病传染模型                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infectedPercent = 0.1;%初始化时网络中存在的感染比例
infectedNumber = floor(userNumber*infectedPercent);
randomNumber = randperm(userNumber);
infectedUsers = randomNumber(1:infectedNumber); %随机选取infectedNumber个节点来作为初始化感染节点
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
%                        选择节点进行免疫                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%1 着色算法锁定用户范围%%%%%%
C = 1:50; %初始化设置20种颜色
colorUsed =  zeros(userNumber,1); % 每个节点使用的颜色向量
saturationDegree = ones(userNumber,1);%节点饱和度
noColorNodes = 1:userNumber; %未染色节点集合
while (~isempty(noColorNodes))
    kk = length(noColorNodes)
    maxSaturationNode = find(saturationDegree == max(saturationDegree(noColorNodes)));
    if length(maxSaturationNode) > 1
        nodeDegreeSet = nodeDegree(maxSaturationNode); %找到最大饱和度的节点集合的度
        colorNode = find(nodeDegree == max(nodeDegreeSet)); %找到最大度的节点
        colorNode = colorNode(1); %找到着色节点
    else
        colorNode = maxSaturationNode;  %找到着色节点
    end
    for c = 1:length(C)
        for j = 1:userNumber
            if A(colorNode,j) > 0 %两个节点相邻
                if colorUsed(j) == c %%节点j使用了颜色c
                    break
                end
            end
        end
        if A(colorNode,j)*colorUsed(j) ~= c %节点j没有用过颜色c
            colorUsed(colorNode) = c;
            break
        end
    end
    noColorNodes(find(noColorNodes==colorNode)) = [];
    saturationDegree(colorNode) = 0;
    nodeDegree(colorNode) = -1;
    for i = 1:userNumber %重新计算饱和度
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
candidateSet = find(colorUsed == maxColor); %%通过着色定理选出了节点
%算法2计算节点重要性：选择容易被传染且容易传染给别人的节点
%1)计算节点的感染概率
%求最短路径，边的最短路径是上面权值的倒数
W_ = W;
for i = 1:userNumber
    for j = (i+1):userNumber
        if W_(i,j) > 0 
            W_(i,j) = 1/W_(i,j);
            W_(j,i) = W_(i,j);
        end
    end
end
infectedScore = zeros(length(candidateSet),1)%candidate的被感染值
%对每个传染源进行处理
for i = 1:length(infectedUsers) %传染源
     i
    mark =DJshortpath(W_,infectedUsers(i)); %求传染源到candidate的最短路径
    for j = 1:length(candidateSet) %%candidate中的节点
        j
         [number,Spth] = findPath(mark,infectedUsers(i),candidateSet(j));     
        shortPathNumber = length(Spth); %最短路径数量
        if shortPathNumber == 1 %%最短路径数量为1
            infected = 1; %%传染概率
            for k = 1:(length(Spth{1})-1) %最短路径节点 
                infected = infected * W(Spth{1}(k),Spth{1}(k+1));
            end
            infectedScore(j)  =  infectedScore(j)+infected; 
        elseif shortPathNumber > 1
            for n = 1:shortPathNumber
                infected = 1;
                for k = 1:(length(Spth{n})-1) %最短路径节点
                    infected = infected * W(Spth{n}(k),Spth{n}(k+1));
                end
                 infectedScore(j)  =  infectedScore(j)+infected; 
            end
        end       
    end
end
%2)求感染值
infectingScore = zeros(length(candidateSet),1);
for i = 1:length(candidateSet)
    shortPathHopSet = cell(userNumber-1,1);
     mark =DJshortpath(W_,candidateSet(i)); %求传染源到candidate的最短路径
     for j = 1:userNumber
         dis=mark{j,2};
         if dis == inf
             break
         end
         [number,Spth] = findPath(mark,candidateSet(i),j);
          shortPathNumber = length(Spth); %最短路径数量
          if shortPathNumber == 1 %%最短路径数量为1
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
     immuneNodeSet = randomNumber(1:immuneNumber); %随机选取infectedNumber个节点来作为初始化感染节点 

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
%比较beita的影响
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
     immuneNodeSet = randomNumber(1:immuneNumber); %随机选取infectedNumber个节点来作为初始化感染节点 

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
              





        
        

            

    
