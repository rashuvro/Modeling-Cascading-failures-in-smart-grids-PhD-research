function commTopo = adj2link(commAdj)
m =  size(commAdj, 1);
n =  size(commAdj, 2);
commTopo = [];
k = 1;
for i = 1 : m
        for j = 1: n
                if commAdj(i, j) == 1
                        commTopo(k, :) = [i j];
                        commAdj(j, i) = 0;
                        k = k + 1;
                end
        end
end

end