% This function assigns probabilities of failure to links based on their
% overload amount and position in the grid and chose one failure for this
% step and assigns some time to it
function [FailedIdx, FailedIndices, moreFailures, LinkProb]=...
    FindMultiFailedLinks(TotalPowerLine,Capacity,mpc1,ListOfFailures,alpha)
            moreFailures=0;
            LinkProb=zeros(length(ListOfFailures),1);
            mu=1; % Rate of failure events (It would be better to define a function  based on state to find the mu)
%             beta=1; % minimum probability of failure for the overloaded line
%            FixedFailProb=0.06;
%             alpha=0.1; % threshold
%%%%% Assign some fixed probability of failure for neighbors of previously
%%%%% failed links
% A=find(ListOfFailures>0);
% for i=1:length(A)
%     k=mpc1.branch(A(i),1); % find one end node
%     m=mpc1.branch(A(i),2);
%     for j=1:length(ListOfFailures)
%         if(i~=j)
%             if((mpc1.branch(j,1)==k) ||(mpc1.branch(j,2)==k) || (mpc1.branch(j,1)==m) || (mpc1.branch(j,2)==m))
%                 LinkProb(j)=FixedFailProb;
%             end
%         end
%     end
% end
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
        if (1-alpha)*Capacity(i)<abs(TotalPowerLine(i)) % If it is overloaded
%             WorkingInterval=(1.5*Capacity(i))-(Capacity(i)-0.01*Capacity(i));
%             currentStatus=abs(TotalPowerLine(i))-(Capacity(i)-0.01*Capacity(i));
%             diffRatio=currentStatus/WorkingInterval;
%             diffRatio=diffRatio*(1-beta);
%             diffRatio=beta+diffRatio;
%             LinkProb(i)=1;
            LinkProb(i)=abs(TotalPowerLine(i)/(1-alpha)*Capacity(i));
            
        end
    end
end
LinkProb(187:end)=0;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Chose among the failed links the link with highest probability of
%%%% failure
FailedIndices = find(listOfFailures>0);

FailedIdx=0;
ProbTest=0;
for i=1:length(FailedIndices)
    if(ProbTest<LinkProb(FailedIndices(i)))
        ProbTest=LinkProb(FailedIndices(i));
        FailedIdx=FailedIndices(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign some random time to the failure event based on the time of the
%%%% last event
% t=exprnd(mu);
% Ftime=currentTime+t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end