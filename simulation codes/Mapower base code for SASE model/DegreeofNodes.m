function degree = DegreeofNodes(topo, NumNodes)
        degree=zeros(1, NumNodes);
        for i = 1 : NumNodes
                degree(i) = 0; % degree of the ith node
                for j=1 : size(topo, 1)
                        if ((topo(j,1) == i) || (topo(j, 2) == i))
                                degree(i) = degree(i) + 1;
                        end
                end
        end
end