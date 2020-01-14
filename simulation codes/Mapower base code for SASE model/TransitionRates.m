close all;
clear all;
clc;

Caps=[12 20 45 48 55 100 130 145 200 325 350 375 400 425 450 500 800 850 1000 2200];

NumData=10;
B=[];
for i=1:NumData
    DataName='States';
    DataName=strcat(DataName, num2str(i));
    load(DataName);
    A=States;
    Concat=[B; A];
    B=Concat;
    clear States;
end
States=B;

MaxTrRate=zeros(length(Cap));
MinTrRate=zeros(length(Cap));
BothSameMax=zeros(length(Cap),1);
BothSameMin=zeros(length(Cap),1);
countMax=zeros(length(Cap),1);
countMin=zeros(length(Cap),1);
for i=1:length(Cap)
    for j=i:length(Cap)
        for k=1:length(States(:,1))
            if(States(k,10)==Cap(i))
                countMax(i)=countMax(i)+1;
                if(States(k,8)~=-1)
                    if(States(k+1,10)~=States(k,10))
                        MaxTrRate(i,find(Cap==States(k+1,10)))=MaxTrRate(i,find(Cap==States(k+1,10)))+1;
                    else
                        if(States(k+1,9)==States(k,9))
                            BothSameMax(i)=BothSameMax(i)+1;
                        end
                    end
                end
            end
        end
    end
end
for i=1:length(Cap)
    for j=1:i
        for k=1:length(States(:,1))
            if(States(k,9)==Cap(i))
                countMin(i)=countMin(i)+1;
                if(States(k,8)~=-1)
                    if(States(k+1,9)~=States(k,9))
                        MinTrRate(i,find(Cap==States(k+1,9)))=MaxTrRate(i,find(Cap==States(k+1,9)))+1;
                    else
                        if(States(k+1,10)==States(k,10))
                            BothSameMin(i)=BothSameMin(i)+1;
                        end
                    end
                end
            end
        end
    end
end

figure(1)
for i=1:length(Cap)
    plot(Caps, MaxTrRate(i,:));
    hold on;
end

figure(2)
for i=1:length(Cap)
    plot(Caps, MinTrRate(i,:));
    hold on;
end