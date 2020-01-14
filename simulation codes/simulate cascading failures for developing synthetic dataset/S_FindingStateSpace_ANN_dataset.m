clc;
clear all;
close all;
define_constants;
CaseName='case300';
FakeCapRate = 1; % fake capacity
%% Parameter initialization
TrueCaps=[50 100 200 400 800]; % quantized capacity
%% Initial load on which the capacity of lines will be determined
mpc1 = loadcase(CaseName);
mpc1 = S_ReadyTheCase(mpc1);
%% keep the original values of number of buses and gens and branches before changing
originalNumBus=length(mpc1.bus(:,1));
originalNumGen=length(mpc1.gen(:,1));
originalNumBran=length(mpc1.branch(:,1));
NumBuses = length(mpc1.bus(:,1));
NumBranches = originalNumBran;
%% Seperate the buses with both load and generators into seperate load and generator buses
[mpc1 LoadGenMatch] = S_SeperateGenAndLoad(mpc1);
%% Calculating the total demand and generation capacity of the grid%%%%%
[WhichInitialLoad, Generation, Demand, DemandIndex]=S_FindFullLoadOfGrid(mpc1);
%%
clear mpc1; % clear because
mpc1 = loadcase(CaseName);
mpc1 = S_ReadyTheCase(mpc1);
%% check if any negative load
for i=1:NumBuses
    if (mpc1.bus(i,3)<0) % if the real power of the bus is negative
        mpc1.bus(i,3) = abs(mpc1.bus(i,3));
    end
end
%% Seperate the buses with both load and generators into seperate load and generator buses
[mpc1 LoadGenMatch] = S_SeperateGenAndLoad(mpc1);
%% Find installed capacity of a transmission line and use it as rateA threshold
[Capacity, FlowCap] = S_CapFinder(WhichInitialLoad,mpc1,TrueCaps,originalNumBran);
mpc1.branch(:,6) = FakeCapRate*Capacity;
CountCaps=zeros(1,length(TrueCaps));
for i=1:length(TrueCaps)
    CountCaps(i)=sum(Capacity==TrueCaps(i));
end
%% Reset mpc1, PD: take mpc with separated load and generator
OriginalMPC = mpc1; % this is the MPC with separated load and generator
clear mpc1;
%% initialize parameters
%DGRatioVector = [0.5, .55, 0.6, .65, 0.7, .75, 0.8, .85, 0.9, 0.95, .99];% r
DGRatioVector = [0.5, .55, 0.6, .65, 0.7, .75, 0.8, .85, 0.9,.95,.99];
len_DGRatioVector = length(DGRatioVector);
%DeltaVector = [0.01 .025, 0.05, .075, 0.1, .125, 0.15, .175,  0.2, .225,  0.25, .275, 0.3, .35, .4, .45, .5];
DeltaVector = [0.01 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,0.35,.4,.45,.5];
len_DeltaVector = length(DeltaVector);% e
%NoCoopPercentageVector = [0.01 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,0.45 0.5, 0.6, 0.7, .8, .9]; % theta
NoCoopPercentageVector = [0.01 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,0.35,.4,.45,.5]; % theta
%NoCoopPercentageVector = [0.01 .025, 0.05, .075, 0.1, .125, 0.15, .175,  0.2, .225,  0.25, .275, 0.3,35,.4]; % theta
len_NoCoopPercentageVector = length(NoCoopPercentageVector);
%FixedFailProbvector = [0, 0.01, 0.02, 0.03, 0.04, 0.05]; % fixed probability of failure for neighbors
%len_FixedFailProbvector = length(FixedFailProbvector);

%%
NumIt = 1000;
%% Generate an initial failure table
iniFailNodes = [];
for i=1:NumIt
    %2 or 3 failures
    b=1+ceil(9*rand);
    randomindex=randperm(NumBranches);
    temp=randomindex(1:b);
  %  iniFailNodes = [iniFailNodes;temp];
    IniFidx=randomindex(1:b);
    
    IniFtable{1,i} = IniFidx;
    clear IniFidx
end
%%
mpc1 = OriginalMPC; % this is the MPC with separated load and generator
BranchMatrix=mpc1.branch;
NumBranches=length(BranchMatrix(:,1));
BusMatrix=mpc1.bus;
NumBuses=length(BusMatrix(:,1));
GenMatrix=mpc1.gen;
NumGens=length(GenMatrix(:,1));
StateCounter = 0; %counts the number of possible states
clear States;
limit = 10000;
States = zeros(limit,18);
tic
for s=1:NumIt
    if StateCounter > limit * 0.99
                break;
    end
    s
    %%  Human error probability
    
    % available time
    x = rand();
    if x <=.42
        PSF_time = 10;
    elseif x>.42 && x<= .97
        PSF_time = 1;
    else
        PSF_time = .1;
    end
    
    %stress
    
    x = rand();
    if x <=.03
        PSF_stress = 5;
    elseif x>.03 && x<= .34
        PSF_stress = 2;
    else
        PSF_stress = 1;
    end
    
    %complexity
    x = rand();
    if x <=.11
        PSF_complexity = 5;
    elseif x>.11 && x<= .78
        PSF_complexity = 2;
    elseif x>.78 && x<= .92
        PSF_complexity  = 1;
    else
        PSF_complexity = 0.1;
    end
    
    %experience
    
    x = rand();
    if x <=.11
        PSF_experience = 10;
    elseif x>.11 && x<= .47
        PSF_experience = 1;
    else
        PSF_experience = .5;
    end
    
    
    %procedures
    x = rand();
    if x <=.06
        PSF_procedures = 50;
    elseif x>.06 && x<= .14
        PSF_procedures = 20;
    elseif x>.14 && x<= .20
        PSF_procedures  = 5;
    else
        PSF_procedures = 1;
    end
    
    %ergonomics
    
    x = rand();
    if x <=.36
        PSF_ergonomics = 50;
    elseif x>.36 && x<= .53
        PSF_ergonomics = 10;
    else
        PSF_ergonomics = 1;
    end
    
    PSF_fitness =1;
    
    %work_process
    
    x = rand();
    if x <=.11
        PSF_work_process = 2;
    elseif x>.11 && x<= .80
        PSF_work_process= 1;
    else
        PSF_work_process = .8;
    end
    
    HEP = (0.01*(PSF_complexity*PSF_ergonomics*PSF_experience*PSF_fitness*PSF_procedures*PSF_stress*PSF_time*PSF_work_process))/ ((0.01*((PSF_complexity*PSF_ergonomics*PSF_experience*PSF_fitness*PSF_procedures*PSF_stress*PSF_time*PSF_work_process)-1))+1);
    HEP_org = HEP;
    DGRatio = DGRatioVector(randi(len_DGRatioVector));
    alpha  = DeltaVector(randi(len_DeltaVector));
    alpha_original = alpha;
    NoCoopPercentage = NoCoopPercentageVector(randi(len_NoCoopPercentageVector));
    NoCoopPercentage_original = NoCoopPercentage;
    
    %% effect of HEP 
    
    if HEP > 0.2 
        NoCoopPercentage = round(NoCoopPercentage+HEP*NoCoopPercentage,2);
        alpha = round(alpha+HEP*alpha,2);
        if NoCoopPercentage>0.4
            NoCoopPercentage =0.4;
        end
        if alpha>0.4
            alpha =0.4;
        end
    end
 
    
  %  FixedFailProb = FixedFailProbvector(randi(len_FixedFailProbvector))
    TotalShed = 0;
    ListOfFailures=zeros(1,NumBranches); % List of failures in one scenario of cascade
    mpc1 = OriginalMPC;  % this is the MPC with separated load and generator
    BranchMatrix=mpc1.branch;
    BusMatrix=mpc1.bus;
    GenMatrix=mpc1.gen;
    NumBranches=length(BranchMatrix(:,1));
    NumBuses=length(BusMatrix(:,1));
    NumGens=length(GenMatrix(:,1));
    % adjoint matrix
    AdjMatrix = zeros(NumBuses,NumBuses);
    for j=1:NumBranches
        AdjMatrix(mpc1.branch(j,1),mpc1.branch(j,2))=1;
        AdjMatrix(mpc1.branch(j,2),mpc1.branch(j,1))=1;
    end
    %% Set the load over the grid to according to the DGRatio (use of r)
            Demand = 0;
            for i=1:length(mpc1.bus(:,1))
                if(mpc1.bus(i,2)==1) % if it is a load bus
                    mpc1.bus(i,3)=mpc1.bus(i,3)*DGRatio*WhichInitialLoad;
                    Demand = Demand + mpc1.bus(i,3); % by PD
                    busDemand(i) = mpc1.bus(i,3); % by PD
                end
            end
     %% use of \theta
            LSCooPercent = ones(length(mpc1.bus(:,1)),1);
            % for all buses
            for i=1:length(mpc1.bus(:,1))
                LSCooPercent(i) = LSCooPercent(i)*(1-NoCoopPercentage);
            end
            % for only geneartors: a filter to make sure gens' values in the vector are 0
            for i=1:length(mpc1.bus(:,1))
                if(mpc1.bus(i,2)==2 || mpc1.bus(i,2)==3)
                    LSCooPercent(i)=0;
                end
            end
            %% seperates the controllable and uncontrollable part of load buses
            % this make dispatchable loads with percentage of theta
            % note that this mpc has different number of genertors than
            % original mpc
            
            mpc1 = S_MakePartialLSContGrid(mpc1,LSCooPercent);
            TotalPowerLine = zeros(1,NumBranches);
            TotalGenDem = zeros(1,NumBuses);
            % Start with failure of k-th link
            NewCapM=0; % Maximum Capacity of the failed lines
            CapSum = 0; % Total capacity of the failed lines
            IniFidx=IniFtable{1,s};
            
            %% Shuvro (find maximum and minimum capcity of the initially
            %%failed lines)
            Cmax =Capacity(IniFidx(1));
            Cmin = Capacity(IniFidx(1));                       
            for i=1:length(IniFidx)-1
                if Capacity(IniFidx(i))>Capacity(IniFidx(i+1))
                    Cmax = Capacity(IniFidx(i));
                    Cmin =Capacity(IniFidx(i+1));
                end
            end
            %% Flow capacity of the initially failed lines
            flowcapacity = 0;
            for i= 1:length(IniFidx)
                flowcapacity = flowcapacity+ FlowCap(IniFidx(i));
            end
          %% installed capacity of the failed lines
          
          installedcapacity = 0;
            for i= 1:length(IniFidx)
                installedcapacity = installedcapacity+ Capacity(IniFidx(i));
            end
          
            
            %%
            initial_number_of_line_failures = length(IniFidx);
            for i=1:length(IniFidx)
                l = IniFidx(i); % index of the failed line
                ListOfFailures(l) = 1; % where ever there is 1 means failure
                % remove the line from adjacency matrix
                [AdjMatrix,mpc1] = S_cutLine(AdjMatrix,mpc1,l);
                
                if Capacity(l) < 5000 % 9900 MW lines are not taken into account
                    CapSum = CapSum + Capacity(l);
                end
                % track the maximum capacity
                if NewCapM<Capacity(l)
                    NewCapM=Capacity(l);
                end
            end
           %% Degree and average path length
            G = graph(AdjMatrix);
            D = degree(G);
            Degree = mean(D);
            path = distances(G);
            for i =1:335
                for j=1:335
                    if (path(i,j) == Inf)
                        path(i,j) =0;
                    end
                end
            end
                        
            sum_row = sum(path);
            total = sum(sum_row);
            NumBuses=length(BusMatrix(:,1));
            APL = total/(NumBuses*(NumBuses-1)); 
            %% Whenever a failure or load shed happens we save the state of the grid
            StateCounter=StateCounter+1;
            States(StateCounter,1)=initial_number_of_line_failures; % total line failure
            States(StateCounter,2)=sum(ListOfFailures); % Maximum Capacity of the failed lines
            States(StateCounter,3)=Cmax; % Amount of load shed happend because of the failure in previous step (To Be assigned)
            States(StateCounter,4)=Cmin; % Amount of load shed in comparing to previous step
           % States(StateCounter,5)=Demand; % Initial load(Demand) over the system
            States(StateCounter,5)=NoCoopPercentage_original; % theta
            States(StateCounter,6)=DGRatio; %r
            States(StateCounter,7)=alpha_original; % capacity estimation error
            States(StateCounter,8)=0; % needs to be fillout based on the next failures that may happen to see if this state is steady state or not
            %States(StateCounter,9)=Capacity(l); % capacity of failed ones
            States(StateCounter,10)=Degree; % Maximum  capacity of failed ones
            States(StateCounter,11)=APL; % Index of the failed one
            States(StateCounter,12)=installedcapacity; % Capacity of the failed one
            States(StateCounter,17)=HEP_org;
            States(StateCounter,15)=Demand; % Time of the failure event
           % States(StateCounter,14)=CapSum; % Accumulation of failed capacities
        %%
            moreFailures=1; % Is any failure happened in previous step?
            counter =0;
            while(moreFailures)
                moreFailures=0; % to see we will have more failures or not
                %%%%%%%%%%%%%%% Find the connected components%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                SG = sparse(AdjMatrix);
                [SS, Components]  = graphconncomp(SG,'DIRECTED',false);
                
                SS;
                numC = zeros(1,SS);
                %% shuvro (for tracking whether there was any island
                %%initially or not
                counter = counter +1;
                island=0;
                if counter ==1
                    if SS>1
                        island =1;
                    end
                end
                %% Run the OPF for every island of the power grid due to failure%%%%%%%%%
                [G,P,VB] = S_islandedGrid(mpc1,Components,SS);
                TotalPowerLine=P;
                TotalGenDem=G;
                %% Check the load shed amount %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
                Served=0;
                generation = 0;
                dispatched_served = 0;
                for i=1:NumBuses
                    if(DemandIndex(i)==1)
                        Served=Served+abs(TotalGenDem(i,1))+abs(TotalGenDem(i,2)); % change by PD
                        dispatched_served = dispatched_served + abs(TotalGenDem(i,1));
                    else
                        generation = generation + abs(TotalGenDem(i,1));
                    end
                end
                TotalShed=round(Demand-Served);
                count_steps =1;
                if count_steps ==1
                    States(StateCounter,13)=TotalShed;
                end
                count_steps =count_steps +1;
                States(StateCounter,18)=generation;
                 States(StateCounter,14)= TotalShed;
                States(StateCounter,16)= Demand - TotalShed;
                
                
                %% Check for more failure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %currentTime=States(StateCounter,13);
                [FailedIndex, moreFailures, LinkProb]=...
                    S_FindFailedLink(TotalPowerLine,Capacity,mpc1,ListOfFailures,alpha);
                link_failprobability = LinkProb;
                %neighbor_fail_prob = FixedFailProb
                if(moreFailures==1)
                  %  StateCounter=StateCounter+1;
                    ListOfFailures(FailedIndex)=1;
                    [AdjMatrix, mpc1]=S_cutLine(AdjMatrix,mpc1,FailedIndex);
                end
                
                if moreFailures==1 % If we have failure we need to save a new state
                    A=find(ListOfFailures>0);
                   % States(StateCounter,1)=initial_number_of_line_failures; % This state has only one total failure in the topology
                    States(StateCounter,2)=sum(ListOfFailures); % Total Capacity of failed ones
                   % States(StateCounter,3)=-1; % Amount of load shed happend because of the failure in previous step (To Be assigned)
                   % States(StateCounter,4)=-1; % Amount of load shed in comparing to previous step
                    States(StateCounter,5)=NoCoopPercentage_original; % Initial load(Demand) over the system
                    States(StateCounter,9)=island; % needs to be fillout based on the next failures that may happen to see if this state is steady state or not
                    %States(StateCounter,10)=Degree; % Index of the failed one
                    %States(StateCounter,12)=Capacity(FailedIndex); % Capacity of the failed one
                   % States(StateCounter,13)=initial_number_of_line_failures; % Time of the failure event
                   
                    % States(StateCounter,15)=Demand;
                    % States(StateCounter,17)=HEP;
                    %  min and max
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
                    %States(StateCounter,3)=MinCap; % min capacity of failed ones
                   % States(StateCounter,2)=MaxCap; % max capacity of failed ones
                   % States(StateCounter-1,8)=StateCounter;
                    NotAbsortState = 1;
                else
                    States(StateCounter,8)=round(flowcapacity); % It means previous state was a steady state
                    NotAbsortState = 0;
                end % end of saving the states
            end
            clear mpc1
            toc
end
