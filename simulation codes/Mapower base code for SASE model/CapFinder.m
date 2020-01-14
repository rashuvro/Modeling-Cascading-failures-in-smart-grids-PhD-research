% This function finds the capacity of the lines
function Capacity=CapFinder(Load,mpc1,TrueCaps,FakeCaps,originalNumBran)
%     define_constants;
    alpha=1; % we multiply this with the power flow in lines to obtain the capacity 
    BusMatrix=mpc1.bus;
    GenMatrix=mpc1.gen;
    BranchMatrix=mpc1.branch;

    NumBuses=length(BusMatrix(:,1));
    NumGens=length(GenMatrix(:,1));
    NumBranches=length(BranchMatrix(:,1));
    
    Capacity=(ones(NumBranches,1)).*9900;
    
    mpc1=DispatchPowerGrid(mpc1,Capacity,Load);

    % New Sizes
    BranchMatrix=mpc1.branch;
    BusMatrix=mpc1.bus;
    GenMatrix=mpc1.gen;

    NumBranches=length(BranchMatrix(:,1));
    NumBuses=length(BusMatrix(:,1));
    NumGens=length(GenMatrix(:,1));

    TotalPowerLine=zeros(1,NumBranches); % Power flow in lines
    TotalGenDem=zeros(1,NumBuses); % Power generation and demand

    AdjMatrix=zeros(NumBuses,NumBuses);
    for j=1:NumBranches
        AdjMatrix(mpc1.branch(j,1),mpc1.branch(j,2))=1;
        AdjMatrix(mpc1.branch(j,2),mpc1.branch(j,1))=1;
    end
    %%%%%%%%%%%%%%%% Find the connected components
    SG = sparse(AdjMatrix);
    [SS, Components]  = graphconncomp(SG,'DIRECTED',false);
    numC=zeros(1,SS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G,P]=islandedGrid(mpc1,Components,SS);
    Capacity=abs(P);
    Capacity=Capacity.*alpha;
    %%%%%%%%%%%New Capacities%%%%%%%%%%%%%%%%%%%%%%%%%
%     FakeCaps = [20 40 60 80 120 180 250 400 800];
    for w=1:length(Capacity)
        for g=1:length(FakeCaps)
            if(g==1)
                if(Capacity(w)<=FakeCaps(g))
                    Capacity(w)=FakeCaps(g);
                    break
                end
            else
                    if(Capacity(w)>FakeCaps(g-1) && Capacity(w)<=FakeCaps(g))
                        Capacity(w)=FakeCaps(g);
                        break
                    end
            end
        end
    end
    
    for w=originalNumBran+1:NumBranches
        Capacity(w)=9900;
    end
    
end
