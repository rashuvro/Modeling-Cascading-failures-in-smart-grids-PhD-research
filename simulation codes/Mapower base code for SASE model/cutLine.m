function [AdjMat, mpc]=cutLine(AdjMatrix,mpc1,FailIndex)            
% Theis number is used for generating failures (by multiplying the
% impedence with them)
              LargeNumber1=100000;
              AdjMatrix(mpc1.branch(FailIndex,1),mpc1.branch(FailIndex,2))=0;
              AdjMatrix(mpc1.branch(FailIndex,2),mpc1.branch(FailIndex,1))=0;
              mpc1.branch(FailIndex,6)=mpc1.branch(FailIndex,6)/mpc1.branch(FailIndex,6);
              mpc1.branch(FailIndex,3)=mpc1.branch(FailIndex,3)*LargeNumber1;
              mpc1.branch(FailIndex,4)=mpc1.branch(FailIndex,4)*LargeNumber1;
              mpc=mpc1;
              AdjMat=AdjMatrix;
end