% This function finds the islands in a grid and solve power flow for each 
% island and return the vector of the flow over the lines (P) and the
% vector of generation and demand in (G)
function [G,P,VB]=islandedGrid(mpc1,Components,SS)
        define_constants;
        
        NumBranches=length(mpc1.branch(:,1));
        NumBuses=length(mpc1.bus(:,1));

        TotalPowerLine=zeros(NumBranches,1); % Power flow in lines
        TotalGenDem=zeros(NumBuses,2); % Power generation and demand
        
        IslandTrack=zeros(SS,2); % Number of components in each island
        % the second column is only useful when we have only one bus in the
        % island and we have the index of that line in the second column
    
        for k=1:SS % Number of islands
            count1=0; % Count number of components in each island
            t=0; % This variable is only useful when only one bus is in island
            for i=1:length(Components)
                if (Components(i)==k)
                    count1=count1+1;
                    t=i; % looks for the only bus in the islands with one bus
                end
            end
            IslandTrack(k,1)=count1;
            IslandTrack(k,2)=t;
        end
        
        for k=1:SS
            if(IslandTrack(k,1)~=1) % when size of island is larger than 1 make a case for that island
                    mpc2 = loadcase('caseBase');  % load an empty case to load it with island information
                    count=0; % Number of elements in the island
                    gcount=0; 
                    brcount=0; % Number of branches in the island
                    Greference=0; % reference exists in that island
                    

                    % List of all the branches and if they have been added to
                    % the island they belong to or not
                    BranchTrack=zeros(1,length(mpc1.branch(:,1))); % keep track of branches that has been added to different islands

                    for i=1:length(Components)
                        if (Components(i)==k)
                                if(mpc1.bus(i,2)==3)
                                        Greference=1; % reference exists in that island
                                end

                            count=count+1; % Number of elements in the island
                            mpc2.bus(count,:)=mpc1.bus(i,:); % put that pus in new case
                            mpc2.bus(count,1)=count; % Bus id in new case
                            BusTracker(count,k)=i; % The bus in island k with id count is equal to orginal bus id i
                            %%%%%% Connect the corresponding links to the bus
                            for l=1:length(mpc1.branch(:,1))
                                if(i==mpc1.branch(l,1)) % if the first end of anyline is equal to that bus
                                    if(Components(mpc1.branch(l,2))==k) % if the other end of the line is in the same category as the other end
                                        if(BranchTrack(l)==0) % if that branch has not been added before
                                            brcount=brcount+1;
                                            mpc2.branch(brcount,:)=mpc1.branch(l,:);
                                            BranchTrack(l)=1;
                                            BRTracker(brcount,k)=l;  % Branch in with id bcount in island k is equal to branch with id l in the origina case
                                        end
                                    end
                                end
                                if(i==mpc1.branch(l,2))% if the second end of anyline is equal to that bus
                                    if(Components(mpc1.branch(l,1))==k) % if the other end of the line is in the same category as the other end
                                        if(BranchTrack(l)==0)% if that branch has not been added before
                                            brcount=brcount+1;
                                            mpc2.branch(brcount,:)=mpc1.branch(l,:);
                                            BranchTrack(l)=1;
                                            BRTracker(brcount,k)=l;  % Branch in with id bcount in island k is equal to branch with id l in the origina case
                                        end
                                    end
                                end
                            end
                            % Till this point the bus and the links connected
                            % to that bus are placed in one island
                            if(mpc1.bus(i,2)==2 || mpc1.bus(i,2)==3) % if that bus is a generator

                                gcount=gcount+1; % count the number of generators
                                for j=1:length(mpc1.gen(:,1))
                                    if(mpc1.gen(j,1)==i) % find the generator correspond to bus i
                                        mpc2.gen(gcount,:)=mpc1.gen(j,:);
                                        mpc2.gen(gcount,1)=count;
                                        mpc2.gencost(gcount,:)=mpc1.gencost(j,:);
                                        GenTracker(gcount,k)=i;
                                    end
                                end
                            end

                        end
                    end

                  %%% Assign new index to the end of lines based on the new
                  %%% indexes for the buses
                    for m=1:count % for every bus in this island
                        for n=1:length(mpc2.branch(:,1))
                            if(mpc2.branch(n,1)==BusTracker(m,k))
                                mpc2.branch(n,1)=m;
                            end
                            if(mpc2.branch(n,2)==BusTracker(m,k))
                                mpc2.branch(n,2)=m;
                            end
                        end
                    end
                    
                %%%%% Now the island is ready%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% Check how much generation capacity and load is in this
                %%%%% island
                    GenCaps=0; % Amount of generators in the island
                    FixedLoad=0; % Amount of fixed (uncontrollabel) load/demand in the island
                    DispatLoad=0; % Amount of dispatchable (controllabel) load/demand in the island
                    if (gcount==0)
                        GenCaps=0;
                    else
                        for g=1:gcount
                            if(mpc2.gen(g,9)>0)
                                GenCaps=GenCaps+mpc2.gen(g,9);
                            end
                        end
                    end
                    for b=1:count
                        if(mpc2.bus(b,3)>0)
                            FixedLoad=FixedLoad+mpc2.bus(b,3);
                        end
                    end
                    for g=1:gcount
                        if(mpc2.gen(g,9)==0)
                            DispatLoad=DispatLoad+abs(mpc2.gen(g,10));
                        end
                    end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% Now consider different scenarios that may happen in
                %%%%%%%% the island %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                NoFlowInIsland=0;  %%% If this is 1 means there will be now flow in the isladn
                if(DispatLoad==0 && FixedLoad==0) %% case when there is no load in the island
                    NoFlowInIsland=1;
                    %%%%%% 
                    for br=1:length(mpc2.branch(:,1))
                        TotalPowerLine( BRTracker(br,k))=-inf;
                    end
                    for bb=1:length(mpc2.bus(:,1))
                        TotalGenDem(BusTracker(bb,k),1)=0;
                        TotalGenDem(BusTracker(bb,k),2)=0;
                    end
                    
                end

                if(GenCaps==0)
                    NoFlowInIsland=1;
                    %%%%%% 
                    for br=1:length(mpc2.branch(:,1))
                        TotalPowerLine( BRTracker(br,k))=-inf;
                    end
                    for bb=1:length(mpc2.bus(:,1))
                        if(BusTracker(bb,k)~=0)  %%%% ??? There may be a probelm here!!!
                            TotalGenDem(BusTracker(bb,k),1)=0;
                            TotalGenDem(BusTracker(bb,k),2)=0;
                        end
                    end
                end

                if(FixedLoad>GenCaps) %%% Fixed load is largr than gen cap
                    NoFlowInIsland=1;
                    %%%%%% ?????? Think of a way to reduce fixed load and
                    %%%%%% have flow
                    for br=1:length(mpc2.branch(:,1))
                        TotalPowerLine( BRTracker(br,k))=-inf;
                    end
                    for bb=1:length(mpc2.bus(:,1))
                        if(BusTracker(bb,k)~=0) %%%% ??? There may be a probelm here!!!
                            TotalGenDem(BusTracker(bb,k),1)=0;
                            TotalGenDem(BusTracker(bb,k),2)=0;
                        else
                            disp('SOME THING WRONG &&&&&&&&&&&&&&&&&&&&&&&&&&&');
                            IslandTrack(k,1)
                            BusTracker
                        end
                    end
                end


                if(NoFlowInIsland==0) % there is flow in the island
                    
                    %%%% each island needs to have a refrences bus
                    if(Greference==0) % There is no reference gen in the island
                        for gr=gcount:-1:1 % chose the first generator to be reference
                            if(mpc2.gen(gr,9)>0)
                                GreferenceI=mpc2.gen(gr,1);
                            end
                        end
                        mpc2.bus(GreferenceI,2)=3;
                    end
                    
                    success=0;
                    opt=mpoption('VERBOSE',0,'OUT_ALL',0,'OUT_BUS',0,'OUT_BRANCH',0,'OUT_ALL_LIM',0);
                    [result2, success]=rundcopf(mpc2,opt);
                    WhatToAddIfNotConverge=mpc2.branch(:,6).*(0.1);
                    while(success==0)
                        mpc2.branch(:,6)=mpc2.branch(:,6)+WhatToAddIfNotConverge;
                        opt=mpoption('VERBOSE',0,'OUT_ALL',0,'OUT_BUS',0,'OUT_BRANCH',0,'OUT_ALL_LIM',0);
                        [result2, success]=rundcopf(mpc2,opt);
                    end
                    
                    for bb=1:length(mpc2.bus(:,1))
                        VB(BusTracker(bb,k))=result2.bus(bb,9);
                    end
                    
                    for br=1:length(mpc2.branch(:,1))
                        TotalPowerLine( BRTracker(br,k))=result2.branch(br, PF);
                    end
                    for bb=1:length(mpc2.bus(:,1))
                        if(mpc2.bus(bb,3)>0)
                            TotalGenDem(BusTracker(bb,k),2)=mpc2.bus(bb,3);
                        end
                        if(mpc2.bus(bb,2)>1) % if the bus is generator
                            for gg=1:length(mpc2.gen(:,1))
                                if(mpc2.gen(gg,1)==bb)
                                    TotalGenDem(GenTracker(gg,k),1)=result2.gen(gg,2);
                                end
                            end
                        end
                    end
                end
            else
                TotalGenDem(IslandTrack(k,2),1)=0;
                TotalGenDem(IslandTrack(k,2),2)=0;
                VB(IslandTrack(k,2))=mpc1.bus(IslandTrack(k,2),9);
            end
            
            clear mpc2;
            clear result2;
            clear success;
        end     
         G=TotalGenDem;
         P=TotalPowerLine;
end