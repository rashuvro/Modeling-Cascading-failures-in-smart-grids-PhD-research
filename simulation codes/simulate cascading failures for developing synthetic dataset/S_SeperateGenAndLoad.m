function [m LoadGenMatch]=S_SeperateGenAndLoad(mpc1)
        define_constants;
        % Find the number of each component
        NumGens=length(mpc1.gen(:,1));
        LoadGenMatch=zeros(NumGens,2);
        count=0; % keep track of number of generators
        
        for i=1:length(mpc1.bus(:,1))
            if(mpc1.bus(i,2)==2 || mpc1.bus(i,2)==3)
                count=count+1;
                if(mpc1.bus(i,3)>0)  % if gen has load also
                    addnew=length(mpc1.bus(:,1))+1;  % new number of bus
                    LoadGenMatch(count,1)=i; % this keep track of corresponding new bus to the generator
                    LoadGenMatch(count,2)=addnew;

                    mpc1.bus(addnew,:)=mpc1.bus(i,:);
                    mpc1.bus(addnew,1)=addnew; % index of new bus
                    mpc1.bus(addnew,2)=1; % make it to be load not gen!          
                    mpc1.bus(addnew,8)=mpc1.bus(addnew,8)*0.99; % Volt Mag
                    mpc1.bus(addnew,9)=mpc1.bus(addnew,9)*0.9; % Volt angle
                    % set the load of old bus to be zero to make it a pure
                    % generator
                    mpc1.bus(i,3)=0; % real power demand
                    mpc1.bus(i,4)=0; % reactive power demand
                    % add new branch for this new load bus
                    addnewBranch=length(mpc1.branch)+1;
                    % connect the new load bus directly to the old gen
                    mpc1.branch(addnewBranch,1)=i; % "from" bus number
                    mpc1.branch(addnewBranch,2)=addnew; % to bus number
                    mpc1.branch(addnewBranch,3)=0.000258; % resistance (p.u.)
                    mpc1.branch(addnewBranch,4)=0.00322; % reactance (p.u.)
                    mpc1.branch(addnewBranch,5)=1.23; % charging susceptance (p.u.)
                    mpc1.branch(addnewBranch,6)=9900; % Rating (here capacity)
                    mpc1.branch(addnewBranch,7)=0;
                    mpc1.branch(addnewBranch,8)=0;
                    mpc1.branch(addnewBranch,9)=0;
                    mpc1.branch(addnewBranch,10)=0;
                    mpc1.branch(addnewBranch,11)=1; % branch status
                    mpc1.branch(addnewBranch,12)=-360;
                    mpc1.branch(addnewBranch,13)=360;
                else
                    LoadGenMatch(count,1)=i;
                    LoadGenMatch(count,2)=-1; % no new load bus added for that bus
                end
            end
        end
        m=mpc1;
        clear mpc1;
end