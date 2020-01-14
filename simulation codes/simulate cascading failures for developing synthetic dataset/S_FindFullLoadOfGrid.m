function [fullLoad, G, D, I]=S_FindFullLoadOfGrid(mpc1)
            define_constants;
            Load=1;
            DGRatio=0;
%%%%%%%%%%%%%%%%%%%% Power Grid Status after failure %%%%%%%%%%%%%%%%%%%%%
            BranchMatrix=mpc1.branch;
            NumBranches=length(BranchMatrix(:,1));
            
            Capacity=(ones(NumBranches,1)).*9900;
            mpc1=S_DispatchPowerGrid(mpc1,Capacity,Load);
           
            BusMatrix=mpc1.bus;
            GenMatrix=mpc1.gen;
            NumBuses=length(BusMatrix(:,1));
            NumGens=length(GenMatrix(:,1));

%%%%% Calculating the total demand and generation capacity of the grid%%%%%
        Demand=0;
        DemandIndex=zeros(1,NumGens);
        Generation=0;
        for i=1:NumGens
            if (mpc1.gen(i,10)<0)
                DemandIndex(mpc1.gen(i,1))=1;
                Demand=Demand + abs(mpc1.gen(i,10));
            else
                Generation=Generation+mpc1.gen(i,9);
            end
        end
        G=Generation;
        D=Demand;
        Load=Generation/Demand;
        I=DemandIndex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 fullLoad=Load;
end