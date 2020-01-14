% This function generates a new case with dispatchable loads
function m=DispatchPowerGrid(mpc1,Capacity,Load)       
       
        BusMatrix=mpc1.bus;
        GenMatrix=mpc1.gen;
        BranchMatrix=mpc1.branch;

        NumBuses=length(BusMatrix(:,1));
        NumGens=length(GenMatrix(:,1));
        NumBranches=length(BranchMatrix(:,1));
        %%%%% Make generators pure generator, their load will be zero %%%%%%%%%%
        for i=1:NumBuses
            if(mpc1.bus(i,2)==2 || mpc1.bus(i,2)==3)
                if(mpc1.bus(i,3)~=0 || mpc1.bus(i,4)~=0)
                    mpc1.bus(i,3)=0;
                    mpc1.bus(i,4)=0;
                end
            end
        end
        counter=0;% keeps track of number of dispatchable loads
        %%%%% Convert demands to dispatchable loads %%%%%%%%%%%%%%%%%%%
        for i=1:NumBuses
            if(mpc1.bus(i,2)==1) % if that bus is a demand
                counter=counter+1;
                mpc1.bus(i,2)=2; % convert it to generator (negative generator/dispatchable loads)
                mpc1.gen(NumGens+counter,1)=i; % add its bus id to gen table
                mpc1.gen(NumGens+counter,2)=-1*(mpc1.bus(i,3)); % Pg
%                 mpc1.gen(NumGens+counter,3)=-1*(mpc1.bus(i,4)); % Qg
                mpc1.gen(NumGens+counter,3)=0; % Qg
                mpc1.gen(NumGens+counter,4)=0; % Qmax
                mpc1.gen(NumGens+counter,5)=-1*(mpc1.bus(i,3)); % Qmin
                mpc1.gen(NumGens+counter,6)=1;
                mpc1.gen(NumGens+counter,7)=100;
                mpc1.gen(NumGens+counter,8)=1;
                mpc1.gen(NumGens+counter,9)=0;
                mpc1.gen(NumGens+counter,10)=-1*Load*(mpc1.bus(i,3)); % Actual demand Pmin       
                mpc1.bus(i,3)=0;
                mpc1.bus(i,4)=0;
            end
        end
        for i=1:NumGens % Cost of generators are all equal as following
            mpc1.gencost(i,1)=1;
            mpc1.gencost(i,2)=0;
            mpc1.gencost(i,3)=0;
            mpc1.gencost(i,4)=2;
            mpc1.gencost(i,5)=0;
            mpc1.gencost(i,6)=0;
            mpc1.gencost(i,7)=12;
            mpc1.gencost(i,8)=240;
        end
        for i=NumGens+1:NumGens+counter % Cost of load shed is as following equal for all
            mpc1.gencost(i,1)=1;
            mpc1.gencost(i,2)=0;
            mpc1.gencost(i,3)=0;
            mpc1.gencost(i,4)=2;
            mpc1.gencost(i,5)=-100;
            mpc1.gencost(i,6)=-3000;
%             mpc1.gencost(i,6)=-2200;
            mpc1.gencost(i,7)=0;
            mpc1.gencost(i,8)=0;
        end
    %%% setting the capacity of lines in the mpc to be considered in the
    %%% optimization
        for i=1:NumBranches
            mpc1.branch(i,6)=Capacity(i);
        end
        m=mpc1;
end