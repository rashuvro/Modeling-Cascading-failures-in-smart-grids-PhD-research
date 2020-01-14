function [AdjMat, mpc]=cutMultiLines(AdjMatrix,mpc1,FailIndices)            
% Theis number is used for generating failures (by multiplying the
% impedence with them)
    for i = 1:length(FailIndices)
        LargeNumber1=100000;
        AdjMatrix(mpc1.branch(FailIndices(i),1),mpc1.branch(FailIndices(i),2))=0;
        AdjMatrix(mpc1.branch(FailIndices(i),2),mpc1.branch(FailIndices(i),1))=0;
        mpc1.branch(FailIndices(i),6)=mpc1.branch(FailIndices(i),6)/mpc1.branch(FailIndices(i),6);
        mpc1.branch(FailIndices(i),3)=mpc1.branch(FailIndices(i),3)*LargeNumber1;
        mpc1.branch(FailIndices(i),4)=mpc1.branch(FailIndices(i),4)*LargeNumber1;
    end
    
    mpc = mpc1;
    AdjMat = AdjMatrix;
end