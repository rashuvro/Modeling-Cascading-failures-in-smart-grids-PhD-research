% This function generates a new case with dispatchable loads
% PD: Convert pure bus (non-generator, only demand) as a negative generator 
% and add them in generator matrix and add cost for this new genartors
% After this all the buses will have associated gen (with +ve and -ve power)
% Reason: Dispatchable loads are modeled as negative gen in Matpower and we
% always minimize generation costs in the OPF
% Output: new genarator matrix, gen cost matrix,
% example: IEEE 118 has 163 buses (after load bus and generation bus separation)
% Among 163 buses 109 load buses and 54 generators buses.
% After dispatched it has 163 generators buses where 54 buses have positve power 
% (old generators) and the rest of the buses with negative power
function m=S_DispatchPowerGrid(mpc1,Capacity,Load)    
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
                mpc1.gen(NumGens+counter,2)=-1*Load*(mpc1.bus(i,3)); % Pg
%                 mpc1.gen(NumGens+counter,3)=-1*(mpc1.bus(i,4)); % Qg
                mpc1.gen(NumGens+counter,3)=0; % Qg
                mpc1.gen(NumGens+counter,4)=0; % Qmax
                mpc1.gen(NumGens+counter,5)=-1*(mpc1.bus(i,3)); % Qmin
                mpc1.gen(NumGens+counter,6)=1; % volt mag
                mpc1.gen(NumGens+counter,7)=100; % base MVA
                mpc1.gen(NumGens+counter,8)=1; % generator status
                mpc1.gen(NumGens+counter,9)=0; % maximum real power, it should be zero since it is -ve generator
                mpc1.gen(NumGens+counter,10)=-1*Load*(mpc1.bus(i,3)); % Actual demand Pmin (minimum real power)       
                mpc1.bus(i,3)=0;
                mpc1.bus(i,4)=0;
            end
        end
        %NumGens+counter;
        
        %
        for i=1:NumGens % Cost of generators are all equal as following
            mpc1.gencost(i,1)=1;
            mpc1.gencost(i,2)=0;
            mpc1.gencost(i,3)=0;
            mpc1.gencost(i,4)=2;
            mpc1.gencost(i,5)=0;
            mpc1.gencost(i,6)=0;
            mpc1.gencost(i,7)=10;
            mpc1.gencost(i,8)=300;
        end
        %}
        
        for i=NumGens+1:NumGens+counter % Cost of load shed is as following equal for all
            mpc1.gencost(i,1)=1;
            mpc1.gencost(i,2)=0;
            mpc1.gencost(i,3)=0;
            mpc1.gencost(i,4)=2;
            mpc1.gencost(i,5)=-100;
            mpc1.gencost(i,6)=-3000;
           % mpc1.gencost(i,6)=-2200;
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