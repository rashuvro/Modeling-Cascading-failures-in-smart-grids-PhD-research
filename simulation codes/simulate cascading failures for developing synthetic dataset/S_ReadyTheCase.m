% This function modifies and cleans the cases that do not have appropriate
% indexes (Bus index)
function mpc=S_ReadyTheCase(mpc1)
        define_constants;
        BusMatrix=mpc1.bus;
        GenMatrix=mpc1.gen;
        BranchMatrix=mpc1.branch;

        NumBuses=length(BusMatrix(:,1));
        NumGens=length(GenMatrix(:,1));
        NumBranches=length(BranchMatrix(:,1));
        for i=1:NumBuses
            BusTrack(i)=mpc1.bus(i,1);
            mpc1.bus(i,1)=i;
        end
        for i=1:NumBranches
            mpc1.branch(i,1)=find(BusTrack==mpc1.branch(i,1));
            mpc1.branch(i,2)=find(BusTrack==mpc1.branch(i,2));
        end
        for i=1:NumGens
            mpc1.gen(i,1)=find(BusTrack==mpc1.gen(i,1));
        end
        mpc=mpc1;
end