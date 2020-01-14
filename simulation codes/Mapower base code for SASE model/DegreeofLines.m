function D=DegreeofLines(NumBranches,mpc)

        degree=zeros(1,NumBranches);
        for i=1:NumBranches
            degree1=0; % one end degree of node
            degree2=0; % other end degree of node
            degree(i)=0; % degree of line
            k=mpc.branch(i,1); % find one end node
            for j=1:NumBranches
                if((mpc.branch(j,1)==k) ||(mpc.branch(j,2)==k))
                    degree1=degree1+1;
                end
            end
            k=mpc.branch(i,2); % find the other end node
            for j=1:NumBranches
                if((mpc.branch(j,1)==k) ||(mpc.branch(j,2)==k))
                    degree2=degree2+1;
                end
            end
            degree(i)=(degree1+degree2-2)/2; % In order to exclude the link itself from calculatiion we add -2
        end
        D=degree;
end