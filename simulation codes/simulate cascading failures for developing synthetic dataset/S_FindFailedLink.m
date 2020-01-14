% This function assigns probabilities of failure to links based on their
% overload amount and position in the grid and chose one failure for this
% step and assigns some time to it
%************* PD made some minor changes ********************************
function [FailedIndex, moreFailures, LinkProb]=...
    S_FindFailedLink(TotalPowerLine,Capacity,mpc1,ListOfFailures,alpha)
moreFailures=0; % initialization
LinkProb=zeros(length(ListOfFailures),1); % initialization
mu=1; % Rate of failure events (It would be better to define a function  based on state to find the mu)
beta=0.0; % minimum probability of failure for the overloaded line
%{
x = randi(5);
if (x ==1)
    FixedFailProbvector = 0.00;
end
if (x ==2)
    FixedFailProbvector  = [0,0,0,0,0,0,0, 0.0, 0.0,0,0, 0.0, 0.0,0,0.0, 0.0, 0.0, 0.0,0, 0.01,0,0,0,0,0,0,0, 0.0, 0.0,0,0, 0.0, 0.0,0,0.0, 0.0, 0.0, 0.0,0, 0.0,0,0,0,0,0,0,0,0,0,0];
end
if (x ==3)
    FixedFailProbvector = [0,0 0.0, 0.0,0, 0.0, 0.0, 0.0, 0.0,0.01];
end
if (x ==4)
    FixedFailProbvector = [0,0, 0.0, 0.0,0,0, 0.0, 0.0,0,0.0, 0.0, 0.0, 0.0,0, 0.01];
end
if (x ==5)
    FixedFailProbvector = [0,0,0,0,0,0,0, 0.0, 0.0,0,0, 0.0, 0.0,0,0.0, 0.0, 0.0, 0.0,0, 0.01];

end
%}
FixedFailProbvector = [0,0,0,0,0,0,0, 0.0, 0.0,0,0, 0.0, 0.0,0,0.0, 0.0, 0.0, 0.0,0, 0.01];
len_FixedFailProbvector = length(FixedFailProbvector);
FixedFailProb = FixedFailProbvector(randi(len_FixedFailProbvector));
%FixedFailProb =0.0001;
%%%%% Assign some fixed probability of failure for neighbors of previously
%%%%% failed links
A=find(ListOfFailures>0);
for i=1:length(A)
    k=mpc1.branch(A(i),1); % find one end node
    m=mpc1.branch(A(i),2);
    for j=1:length(ListOfFailures)
        if(i~=j)
            if((mpc1.branch(j,1)==k) ||(mpc1.branch(j,2)==k) || (mpc1.branch(j,1)==m) || (mpc1.branch(j,2)==m))
                LinkProb(j)=FixedFailProb;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Assign prob of failure of failed lines based on their amount of
%%%%% overload
% for i=1:length(ListOfFailures)
%     diffRatio=0;
%     if(TotalPowerLine(i)~=-inf) % it is not a island without any power
%         if ((Capacity(i)-0.01*Capacity(i))<abs(TotalPowerLine(i))) % If it is overloaded
%             WorkingInterval=(1.5*Capacity(i))-(Capacity(i)-0.01*Capacity(i));
%             currentStatus=abs(TotalPowerLine(i))-(Capacity(i)-0.01*Capacity(i));
%             diffRatio=currentStatus/WorkingInterval;
%             diffRatio=diffRatio*(1-beta);
%             diffRatio=beta+diffRatio;
%             LinkProb(i)=LinkProb(i)+diffRatio;
%         end
%     end
% end
for i=1:length(ListOfFailures)
    %     diffRatio=0;
    if(TotalPowerLine(i)~=-inf) % it is not a island without any power
        %if (1-alpha)*Capacity(i)<round(abs(TotalPowerLine(i))) % Changed by PD
        if (1-alpha)*(Capacity(i)+Capacity(i)*0.2)<abs(TotalPowerLine(i)) % If it is overloaded    
            %             WorkingInterval=(1.5*Capacity(i))-(Capacity(i)-0.01*Capacity(i));
            %             currentStatus=abs(TotalPowerLine(i))-(Capacity(i)-0.01*Capacity(i));
            %             diffRatio=currentStatus/WorkingInterval;
            %             diffRatio=diffRatio*(1-beta);
            %             diffRatio=beta+diffRatio;
            %             LinkProb(i)= LinkProb(i)+diffRatio;
            % LinkProb(i)=((1-alpha)*Capacity(i))/ abs(TotalPowerLine(i));
            % ************* by PD from Mahshid's code ********************
            LinkProb(i) = abs(TotalPowerLine(i))/Capacity(i); 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Find the failed links based on the probabilities obtained above
listOfFailures=zeros(length(ListOfFailures),1);
for i=1:length(ListOfFailures)
    if(rand<LinkProb(i))
        if(ListOfFailures(i)~=1)
            listOfFailures(i)=1;
            moreFailures=1;
        end
    end
end

%% ************* changed by PD from intuition ********************
%{
for i=1:length(ListOfFailures)
    if(LinkProb(i)>0)
        if(ListOfFailures(i)~=1)
            listOfFailures(i)=1;
            moreFailures=1;
        end
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Chose among the failed links the link with highest probability of
%%%% failure

B=find(listOfFailures>0);
%{
if length(B) == 0
    FailedIndex=0;
else
    x =randi(length(B));
    FailedIndex =B(x);
end
%}

FailedIndex=0;
ProbTest=0;
for i=1:length(B)
    if(ProbTest<LinkProb(B(i)))
        ProbTest=LinkProb(B(i));
        FailedIndex=B(i);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign some random time to the failure event based on the time of the
%%%% last event
% t=exprnd(mu);
% Ftime=currentTime+t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end