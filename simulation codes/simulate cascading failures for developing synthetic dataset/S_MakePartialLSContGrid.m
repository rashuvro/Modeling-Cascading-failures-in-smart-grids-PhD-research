% seperates the controllable and uncontrollable part of load buses
% controllable parts are converted to generators with negative generations
% uncontrollable parts are remained in the bus as the fixed load
% FYI, this function has almost similar operation as dispatched load
function m=S_MakePartialLSContGrid(mpc1,LSCVector)

NumGens=length(mpc1.gen(:,1));

for i=1:length(mpc1.bus(:,1))
    if(LSCVector(i)~=0) % if this is not a pure generator (this is zero for pure generators)
        AddNewGen=length(mpc1.gen(:,1))+1;
        mpc1.bus(i,2)=2; % convert it to generator (negative generator/dispatchable loads)
        mpc1.gen(AddNewGen,1)=i; % add its bus id to gen table
        mpc1.gen(AddNewGen,2)= -1*LSCVector(i)*(mpc1.bus(i,3)); % Pg (Assign the load on the bus as negative gen value)
        mpc1.gen(AddNewGen,3)= - 1*(mpc1.bus(i,4)); % Qg
        mpc1.gen(AddNewGen,4)= 0; % Qmax
        mpc1.gen(AddNewGen,5)= - 1*min(abs(mpc1.bus(i,4)),1);
        %mpc1.gen(NumGens+counter,5)=-1*(mpc1.bus(i,3)); % Qmin, we assign it to be same as real Pg
        mpc1.gen(AddNewGen,6)=1;
        mpc1.gen(AddNewGen,7)=100;
        mpc1.gen(AddNewGen,8)=1;
        mpc1.gen(AddNewGen,9)=0;
        mpc1.gen(AddNewGen,10)= -1*LSCVector(i)*(mpc1.bus(i,3)); % Actual demand Pmin
        mpc1.bus(i,3)= (1-LSCVector(i))*(mpc1.bus(i,3)); % this was not done in dispacthed load function
        mpc1.bus(i,4)= 0;
    end
end

% this cost is not included in Mahshid' code
 %{
for i=1:NumGens % Cost of generators are all equal as following
    mpc1.gencost(i,1)=1;
    mpc1.gencost(i,2)=0;
    mpc1.gencost(i,3)=0;
    mpc1.gencost(i,4)=2;
    mpc1.gencost(i,5)=0;
    mpc1.gencost(i,6)=0;
    mpc1.gencost(i,7)=10;
    mpc1.gencost(i,8)=100;
end
%}
for i=NumGens+1:length(mpc1.gen(:,1)) % Cost of load shed is as following equal for all
    mpc1.gencost(i,1)=1;
    mpc1.gencost(i,2)=0;
    mpc1.gencost(i,3)=0;
    mpc1.gencost(i,4)=2;
    mpc1.gencost(i,5)=-100;
    mpc1.gencost(i,6)=-4300;
    % mpc1.gencost(i,6)=-2200;
    mpc1.gencost(i,7)=0;
    mpc1.gencost(i,8)=0;
end
m=mpc1;
clear mpc1;
end