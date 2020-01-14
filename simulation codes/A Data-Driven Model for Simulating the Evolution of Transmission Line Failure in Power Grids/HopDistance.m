function H=HopDistance(NumBranches,NumBuses,mpc)
        distance=zeros(NumBranches,NumBranches); % Hop distance of every two links
        AdjMatrix=zeros(NumBuses,NumBuses);
            for j=1:NumBranches
                AdjMatrix(mpc.branch(j,1),mpc.branch(j,2))=1;
                AdjMatrix(mpc.branch(j,2),mpc.branch(j,1))=1;
            end
            for i=1:NumBuses
                for j=1:NumBuses
                    if (j==i)
                         WeiMatrix(i,i)=0;
                    else
                        if(AdjMatrix(i,j)==1)
                            WeiMatrix(i,j)=1;
                        end
                        if(AdjMatrix(i,j)==0)
                            WeiMatrix(i,j)=inf;
                        end
                    end
                end
            end
            AllSpSf = sp(AdjMatrix,WeiMatrix); % Distance of every pair of nodes in the system
            
        for i=1:NumBranches
           for j=1:NumBranches
                  d1=AllSpSf(mpc.branch(i,1),mpc.branch(j,1));
                  d2=AllSpSf(mpc.branch(i,2),mpc.branch(j,1));
                  d3=AllSpSf(mpc.branch(i,1),mpc.branch(j,2));
                  d4=AllSpSf(mpc.branch(i,2),mpc.branch(j,2));
                  distance(i,j)=min(min(min(d1,d2),d3),d4);
           end
        end
        H=distance;
end