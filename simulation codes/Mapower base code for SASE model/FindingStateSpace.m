clc;
clear all;
close all;
define_constants;
CaseName='case118';
FakeCapRate =1;
% Parameter initialization
TrueCaps=[20 80 200 500 800]; % quantized capacity 
FakeCaps = TrueCaps;
NumIt = 10; % Number of iteration for extracting states..
% Initial load on which the capacity of lines will be determined
mpc1 = loadcase(CaseName);
mpc1=ReadyTheCase(mpc1);
% Seperate the buses with both load and generators into seperate load and generator buses 
[mpc1 LoadGenMatch]=SeperateGenAndLoad(mpc1);
[WhichInitialLoad, Generation, Demand, DemandIndex]=FindFullLoadOfGrid(mpc1); % In the case of IEEE 118 topology 6.9 corresponds to using the full generation capacity of the grid
%WhichInitialLoad = 0.2537
clear mpc1; % clear because 
mpc1 = loadcase(CaseName);
mpc1=ReadyTheCase(mpc1);
% keep the original values of number of buses and gens and branches before changing
originalNumBus=length(mpc1.bus(:,1));
originalNumGen=length(mpc1.gen(:,1));
originalNumBran=length(mpc1.branch(:,1));
NumBuses = length(mpc1.bus(:,1));
% check if any negative load
for i=1:NumBuses
    if (mpc1.bus(i,3)<0)
        mpc1.bus(i,3) = abs(mpc1.bus(i,3));
    end
end
% Seperate the buses with both load and generators into seperate load and generator buses
[mpc1 LoadGenMatch]=SeperateGenAndLoad(mpc1);
% Find installed capacity of a transmission line and use it as rateA threshold
Capacity=CapFinder(WhichInitialLoad,mpc1,TrueCaps,FakeCaps,originalNumBran);

mpc1.branch(:,6)=Capacity;
CountCaps=zeros(1,length(TrueCaps));
for i=1:length(TrueCaps)
    CountCaps(i)=sum(Capacity==TrueCaps(i));
end
tic
% Reset mpc1
OriginalMPC=mpc1;
NumBranches=231;
% Generate a initial failure table
for i=1:NumIt
   %{ % 2 or 3 failures
  
 b=1+ceil(2*rand);
     b=3;
    randomindex=randperm(NumBranches);
    temp=randomindex(1:b);
    IniFidx=randomindex(1:b);
   IniFtable{1,i} = IniFidx;
   clear IniFidx
end
NoCoopPercentageVector = 0.1;
%DGRatioVector=[0.01 0.5 0.9 0.99];
DGRatioVector = 0.9;
DeltaVector = 0.1;
% DeltaVector=[0.1 0.2 0.3 0.4 0.45 0.5];

tic
for idx=1:length(NoCoopPercentageVector)
    NoCoopPercentage=NoCoopPercentageVector(idx);
    for DGRatio=DGRatioVector
        clear mpc1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mpc1 = OriginalMPC;
            BranchMatrix=mpc1.branch;
            NumBranches=length(BranchMatrix(:,1));
            BusMatrix=mpc1.bus;
            NumBuses=length(BusMatrix(:,1));
            GenMatrix=mpc1.gen;
            NumGens=length(GenMatrix(:,1));

            Degrees = DegreeofLines(NumBranches,mpc1); % Find the degree of line
            HopDist = HopDistance(NumBranches,NumBuses,mpc1);% Finds the hop distance between all links
            
            %%%%% Calculating the total demand and generation capacity of the grid%%%%%
            %         =Demand/Generation;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        for alpha=DeltaVector
                clear States;
                StateCounter=0; %counts the number of possible states
                
                limit = 10000;
                States = zeros(limit, 13);

                for s=1:NumIt % for every iteration
                    if StateCounter > limit * 0.99
                        break;
                    end
                    ListOfFailures=zeros(1,NumBranches); % List of failures in one scenario of cascade
                    
                %%%%%%%%%%%%% convert the loads in the power grid to dispatchable loads
                    mpc1 = OriginalMPC;

                    for i=1:length(mpc1.bus(:,1))
                        if(mpc1.bus(i,2)==1)
                            mpc1.bus(i,3)=mpc1.bus(i,3)*DGRatio*WhichInitialLoad;
                        end
                    end
                 %%%%%%%%%%%%%%%%%%%% Power Grid Status after failure %%%%%%%%%%%%%%%%%%%%%
                    BranchMatrix=mpc1.branch;
                    BusMatrix=mpc1.bus;
                    GenMatrix=mpc1.gen;
                    NumBranches=length(BranchMatrix(:,1));
                    NumBuses=length(BusMatrix(:,1));
                    NumGens=length(GenMatrix(:,1));


                    LSCooPercent=ones(length(mpc1.bus(:,1)),1);
                    for i=1:length(mpc1.bus(:,1))
                        LSCooPercent(i)=LSCooPercent(i)*(1-NoCoopPercentage);
                    end
                    %%% A filter to make sure gens' values in the vector are 0
                    for i=1:length(mpc1.bus(:,1))
                        if(mpc1.bus(i,2)==2 || mpc1.bus(i,2)==3)
                            LSCooPercent(i)=0;
                        end
                    end

                    mpc1=MakePartialLSContGrid(mpc1,LSCooPercent);
                    mpc1.branch(:,6)=2*Capacity;
                    TotalPowerLine=zeros(1,NumBranches);
                    TotalGenDem=zeros(1,NumBuses);

                    AdjMatrix=zeros(NumBuses,NumBuses);
                    for j=1:NumBranches
                        AdjMatrix(mpc1.branch(j,1),mpc1.branch(j,2))=1;
                        AdjMatrix(mpc1.branch(j,2),mpc1.branch(j,1))=1;
                    end

                    %%%%%%%%%%%%%%%%%% Initial failure, Add failures %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    IniFidx=IniFtable{1,s};
                    % Start with failure of k-th link
                    NewCapM=0;
                    CapSum = 0;
                    for i=1:length(IniFidx)
                        k=IniFidx(i);
                        if Capacity(k) < 5000
                            CapSum = CapSum + Capacity(k);
                        end
                        [AdjMatrix, mpc1]=cutLine(AdjMatrix,mpc1,k);
                        ListOfFailures(k)=1; % where ever there is 1 means failure
                        if NewCapM<Capacity(k)
                            NewCapM=Capacity(k);
                        end
                    end
                    
                      %Whenever a failure or load shed happens we save the
                      %state of the grid
                      StateCounter=StateCounter+1;
                      States(StateCounter,1)=sum(ListOfFailures); % This state has only one total failure in the topology
                      States(StateCounter,2)=NewCapM; % Total Capacity of failed ones
                      States(StateCounter,3)=-1; % Amount of load shed happend because of the failure in previous step (To Be assigned)
                      States(StateCounter,4)=-1; % Amount of load shed in comparing to previous step
                      States(StateCounter,5)=Demand; % Initial load(Demand) over the system
                      States(StateCounter,6)=Degrees(k); % degree of the failed link
                      States(StateCounter,7)=HopDist(k,k); % Average hop distance of failures
                      States(StateCounter,8)=0; % needs to be fillout based on the next failures that may happen to see if this state is steady state or not
                      States(StateCounter,9)=Capacity(k); % min capacity of failed ones
                      States(StateCounter,10)=Capacity(k); % max capacity of failed ones
                      States(StateCounter,11)=k; % Index of the failed one
                      States(StateCounter,12)=Capacity(k); % Capacity of the failed one
                      States(StateCounter,13)=0; % Time of the failure event
                      States(StateCounter,14)=CapSum; % Accumulation of failed capacities
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                moreFailures=1; % Is any failure happened in previous step?
                    while(moreFailures)
                        moreFailures=0; % to see we will have more failures or not
                %%%%%%%%%%%%%%%% Find the connected components%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        SG = sparse(AdjMatrix);
                        [SS, Components]  = graphconncomp(SG,'DIRECTED',false);
                        numC=zeros(1,SS);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% Run the OPF for every island of the power grid due to failure%%%%%%%%%
                        [G,P,VB]=islandedGrid(mpc1,Components,SS);
                        TotalPowerLine=P;
                        TotalGenDem=G;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%% Check the load shed amount %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Served=0;
                        for i=1:NumBuses
                            if(DemandIndex(i)==1)
                                Served=Served+abs(TotalGenDem(i));
                            end
                        end
                        TotalShed=round(Demand-Served);
                        States(StateCounter,3)=TotalShed;
                        if(States(StateCounter,13)==0) % if this is the first failure
                            States(StateCounter,4)=TotalShed;
                        else
                            States(StateCounter,4)=States(StateCounter,3)-States(StateCounter-1,3);
                        end
                %%%%% Set the new generation values as the initial gen and demand of buses%
                    for i=1:NumGens
                        mpc1.gen(i,2)=TotalGenDem(i);
                    end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%% Check for more failure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        currentTime=States(StateCounter,13);
                        [FailedIndex, moreFailures, LinkProb]=...
                        FindFailedLink(TotalPowerLine,2*Capacity,mpc1,ListOfFailures,alpha);
                        if(moreFailures==1)
                            StateCounter=StateCounter+1;
                            ListOfFailures(FailedIndex)=1;
                            [AdjMatrix, mpc1]=cutLine(AdjMatrix,mpc1,FailedIndex);
                        end

                        if moreFailures==1 % If we have failure we need to save a new state
                                A=find(ListOfFailures>0);
                                Cap=0;
                                Deg=0;
                                Hop=0;
                                for m=1:length(A)
                                    Cap=Cap+Capacity(A(m));
                                    Deg=Deg+Degrees(A(m))/length(A);
                                    for n=1:length(A)
                                        Hop=Hop+HopDist(m,n)/length(A);
                                    end
                                end
                              States(StateCounter,1)=sum(ListOfFailures); % This state has only one total failure in the topology
                              States(StateCounter,2)=Cap; % Total Capacity of failed ones
                              States(StateCounter,3)=-1; % Amount of load shed happend because of the failure in previous step (To Be assigned)
                              States(StateCounter,4)=-1; % Amount of load shed in comparing to previous step
                              States(StateCounter,5)=Demand; % Initial load(Demand) over the system
                              States(StateCounter,6)=Deg; % degree of the failed link
                              States(StateCounter,7)=Hop; % Average hop distance of failures
                              States(StateCounter,8)=0; % needs to be fillout based on the next failures that may happen to see if this state is steady state or not

                              States(StateCounter,11)=FailedIndex; % Index of the failed one
                              States(StateCounter,12)=Capacity(FailedIndex); % Capacity of the failed one
    %                           States(StateCounter,13)=Ftime; % Time of the failure event
                              States(StateCounter,14)=States(StateCounter-1,14)+States(StateCounter,12); % Accumulation of failed capacities
                                %min and max
                                MinCap=max(TrueCaps);
                                MaxCap=0;
                                for g=1:NumBranches
                                    if (ListOfFailures(g)==1)
                                        if(Capacity(g)<=MinCap)
                                            MinCap=Capacity(g);
                                        end
                                        if(Capacity(g)>=MaxCap)
                                            MaxCap=Capacity(g);
                                        end
                                    end
                                end
                                States(StateCounter,9)=MinCap; % min capacity of failed ones
                                States(StateCounter,10)=MaxCap; % max capacity of failed ones
                                States(StateCounter-1,8)=StateCounter;
                        else
                            States(StateCounter,8)=-1; % It means previous state was a steady state
                        end
                    end
                end
        %         save(strcat('States2012-OPFloading',num2str(100*DGRatio),'.mat'),'States');
%                save(strcat('States2016-OPFother',num2str(1000*DGRatio),'Alpha',num2str(alpha),...
%                   'Beta',num2str(NoCoopPercentage), '.mat'), 'States');
            %     save(strcat('States2012-OPF07Beta',num2str(beta),'.mat'),'States');
        end
    end
end
toc
Blackoutsize=zeros(NumBranches,1);
for i=1:length(States(:,1))
    if(States(i,8)==-1)
        Blackoutsize(States(i,1))=Blackoutsize(States(i,1))+1;
    end
end
