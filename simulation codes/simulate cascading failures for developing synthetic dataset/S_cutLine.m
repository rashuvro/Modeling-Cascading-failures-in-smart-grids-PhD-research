% this function removes the transmission line given by FailIndex from
% Adjacency matrix and MPC
function [AdjMat, mpc]=S_cutLine(AdjMatrix,mpc1,FailIndex)

LargeNumber1=100000; % This number is used for generating failures 
                     % by multiplying the impedence with them)
% change adjance matrix
AdjMatrix(mpc1.branch(FailIndex,1),mpc1.branch(FailIndex,2))=0;
AdjMatrix(mpc1.branch(FailIndex,2),mpc1.branch(FailIndex,1))=0;
% change the MPC
mpc1.branch(FailIndex,3)=mpc1.branch(FailIndex,3)*LargeNumber1; % resistance
mpc1.branch(FailIndex,4)=mpc1.branch(FailIndex,4)*LargeNumber1; % reactance
% line capacity rating
mpc1.branch(FailIndex,6)=mpc1.branch(FailIndex,6)/mpc1.branch(FailIndex,6);
mpc=mpc1;
AdjMat=AdjMatrix;
end