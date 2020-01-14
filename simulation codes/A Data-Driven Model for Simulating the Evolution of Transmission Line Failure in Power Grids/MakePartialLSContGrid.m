% seperates the controllable and uncontrollable part of load buses
% it makes them as a generator with some fixed load
function m=MakePartialLSContGrid(mpc1, LSCVector)
        
        NumGens=length(mpc1.gen(:,1));
        
        for i=1:length(mpc1.bus(:,1))
            if(LSCVector(i)~=0)
                AddNewGen=length(mpc1.gen(:,1))+1;
                mpc1.bus(i,2)=2; % convert it to generator (negative generator/dispatchable loads)
                mpc1.gen(AddNewGen,1)=i; % add its bus id to gen table
                mpc1.gen(AddNewGen,2)=-1*0.5*LSCVector(i)*(mpc1.bus(i,3)); % Pg (Assign the load on the bus as negative gen value)
                mpc1.gen(AddNewGen,3)=-1*(mpc1.bus(i,4)); % Qg
                mpc1.gen(AddNewGen,4)=0; % Qmax
                mpc1.gen(AddNewGen,5)=-1*min(abs(mpc1.bus(i,4)),1);
                %mpc1.gen(NumGens+counter,5)=-1*(mpc1.bus(i,3)); % Qmin, we assign it to be same as real Pg
                mpc1.gen(AddNewGen,6)=1;
                mpc1.gen(AddNewGen,7)=100;
                mpc1.gen(AddNewGen,8)=1;
                mpc1.gen(AddNewGen,9)=0;
                mpc1.gen(AddNewGen,10)=-1*LSCVector(i)*(mpc1.bus(i,3)); % Actual demand Pmin       
                mpc1.bus(i,3)=(1-LSCVector(i))*(mpc1.bus(i,3));
                mpc1.bus(i,4)=0;
            end
        end
        
        % we assign equal cost of load shedding for all loads
        for i=NumGens+1:length(mpc1.gen(:,1)) % Cost of load shed is as following equal for all
            mpc1.gencost(i,1)=1;
            mpc1.gencost(i,2)=0;
            mpc1.gencost(i,3)=0;
            mpc1.gencost(i,4)=2;
            mpc1.gencost(i,5)=-100;
            mpc1.gencost(i,6)=-3000;
            mpc1.gencost(i,7)=0;
            mpc1.gencost(i,8)=0;
        end
        m=mpc1;
        clear mpc1;
end