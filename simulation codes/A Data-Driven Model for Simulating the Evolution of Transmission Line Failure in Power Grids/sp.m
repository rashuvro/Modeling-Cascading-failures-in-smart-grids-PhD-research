%FUNCTION TO COMPUTE THE SHORTEST PATH 
%Input arguments:
%   A:  Adjacency Matrix
%   E:  Matrix of Euclidean distances 

function AllSpSf2 = sp(A,E)
NumNodes2=length(A);
Weight2=zeros(NumNodes2, NumNodes2);
D2=zeros(NumNodes2,NumNodes2,NumNodes2);
AllSpSf2=zeros(NumNodes2,NumNodes2);
for i=1:NumNodes2
    for j=1:i-1
        if(A(i,j)==1)
            Weight2(i,j)=E(i,j);
        else
            Weight2(i,j)=inf;
        end
    end
    
end
Weight2=Weight2+Weight2.';

D2(1,:,:)=Weight2;
for k=1:NumNodes2
    for i=1:NumNodes2
        for j=1:NumNodes2
            if(k>1)
                D2(k,i,j)=min(D2(k-1,i,j),D2(k-1,i,k)+D2(k-1,k,j));
            else
                D2(k,i,j)=min(Weight2(i,j),Weight2(i,k)+Weight2(k,j));
            end
        end 
    end 
end

for i=1:NumNodes2
    for j=1:NumNodes2
        AllSpSf2(i,j)=D2(NumNodes2,i,j);
    end
end

