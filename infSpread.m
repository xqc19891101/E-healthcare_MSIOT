function dy = infSpread(t,y)
global W
global beita
global mu
userNumber = size(W,1);
dy = zeros(3*userNumber,1);
for i = 1:userNumber
    sum = 0;
    for j = 1:userNumber
        sum = W(i,j) * y(userNumber+j) + sum;
    end
    dy(i) = -beita * y(i) * sum;
    dy(userNumber+i) = beita * y(i) * sum - y(userNumber+i)*mu;
    dy(2*userNumber+i) = y(userNumber+i)*mu;
end
    
    