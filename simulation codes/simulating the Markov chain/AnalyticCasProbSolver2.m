% Conditional Blackout Size distribution for different destinations
% Input parameter table
clear all; close all;clc;
load tail.mat
% IEEE-118 bus system
NumberOfLines=186;
% Truncation (trick to reduce numerical solution time)
L=186;
epsilon=0.01;
% Max number of capacities
Capa=[20 80 200 500 800];
% Caps=[80 200 500 800 1500];
C=length(Capa);
% Choose initial state INDEX
Fi=2; Ci=1;

ii=2*C*(Fi-1)+2*Ci-1;

% Weights
Wf=0.5; Wcmax=0.5;
deltaT=0.1;

tic
ParameterTable2V2
Dist2=zeros(length(Level),length(Alpha),NumberOfLines);
for m=1:length(Level)
    for n=1:length(Alpha)
        % Preallocate transition matrix with size 2*NumberOfLines*C
        Q=zeros(2*NumberOfLines*C,2*NumberOfLines*C);
    %     parameter=PARAMETER(m);
    %     IndexParameterReplace
%         a1=A1(m,n);
%         a2=A2(m,n);
%         a3=A3(m,n);
%         a4=A4(m,n);
        a1=ParaSetting2{m,n}(1);
        a2=ParaSetting2{m,n}(2);
        a3=ParaSetting2{m,n}(3);
        a4=ParaSetting2{m,n}(4);
        % If start from 1 failure
        for i=1:NumberOfLines-1
            % j==1 for transient, j==2 for stable states
            % find stability probability parametricly first
            % P_{stop}(F_i)
            if i<=floor(a2*NumberOfLines)
                f1=epsilon + a1*( (a2*NumberOfLines-i)/(a2*NumberOfLines) )^4;
            end
            if i>floor(a2*NumberOfLines) && i<=floor(0.6*NumberOfLines)
                f1=epsilon;
            end
            if i>floor(0.6*NumberOfLines)
                f1=min(1, (epsilon + ( (i-0.6*NumberOfLines)/...
                    (NumberOfLines-0.6*NumberOfLines) )^4) );
            end
            for j=1:2
                % For transiant states j==1
                if j==1
                    % If failure starts from smallest capacity
                    for k=1:C
                        % P_{stop}(C^{\max}_i)
                        f2=max(a4, a3*( (Capa(k)-max(Capa))/max(Capa) )^4 );
                        % Weighted average
                        pStable=Wf*f1 + Wcmax*f2;
                        if pStable>1
                            pStable=1;
                        end
%                         pStableM(i,k)=pStable;
                        % Continue probability
                        pCont=1-pStable;
                        % Rates for cascade-stop transitions
                        Q( 2*(i-1)*C+2*(k-1)+j, 2*(i-1)*C+2*(k-1)+j+1)=pStable;
                        % Rates for cascade-continue transitions
                        % Going to same Cmax
                        if k==1
                            pContSub=min(1,0.03+6e-7*(i+112)^3);
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+1)=...
                                pCont*(1-pContSub);
                        end
                        if k==2
                            pContSub=min(1,0.03+6e-7*(i+75)^3);
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+1)=...
                                pCont*(1-pContSub);
                        end
                        if k==3
                            pContSub=min(1,0.03+6e-7*(i+20)^3);
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+1)=...
                                pCont*(1-pContSub);
                        end
                        if k==4
                            pContSub=min(1,0.03+6e-7*(i-60)^3);
                            if pContSub<0.03
                                pContSub=0.03;
                            end
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+1)=...
                                pCont*(1-pContSub);
                        end
                        if k==5
                            pContSub=0;
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+1)=...
                                pCont*(1-pContSub);
                        end
                        
                        % Going to larger Cmax
                        a=2.22;
                        b=1.52;
                        c=0.52;
                        d=0.03;
                        w=[a b c d];
                        if k==1
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+3)=...
                                pContSub*pCont*w(1)/(sum(w));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+5)=...
                                pContSub*pCont*w(2)/(sum(w));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+7)=...
                                pContSub*pCont*w(3)/(sum(w));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+9)=...
                                pContSub*pCont*w(4)/(sum(w));
                        end
                        if k==2
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+3)=...
                                pContSub*pCont*w(2)/(sum(w)-w(1));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+5)=...
                                pContSub*pCont*w(3)/(sum(w)-w(1));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+7)=...
                                pContSub*pCont*w(4)/(sum(w)-w(1));
                        end
                        if k==3
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+3)=...
                                pContSub*pCont*w(3)/(sum(w)-w(1)-w(2));
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+5)=...
                                pContSub*pCont*w(4)/(sum(w)-w(1)-w(2));
                        end
                        if k==4
                            Q(2*(i-1)*C+2*(k-1)+1,2*i*C+2*(k-1)+3)=...
                                pContSub*pCont*1;
                        end
                    end
                end
            end
        end
        % Assign rates=0 for F_i=NumberOfLines
        Q(2*(NumberOfLines-1)*C+1:2*C*NumberOfLines,:)=0;
        for i=1:2*NumberOfLines*C
            Q(i,i)=-1*sum(Q(i,:));
        end
        Q=Q/deltaT;

%         Q1=Q(1:L*C*2,1:L*C*2);
%         Q1(end-39:end,end-39:end)=qt;

        [V,D] = eig(Q);
        U=diag(V);
        
        M=1;
        for j=2:2:2*NumberOfLines*C
            % Analytical
            if U(j)~=0
                temp=V(ii,j)/U(j);
                Dist2(m,n,M)=Dist2(m,n,M)+temp;
            end
            if mod(j,2*C)==0
                M=M+1;
            end
        end
    end
end
toc
save('Dist2V2.mat','Dist2','Level','Alpha','-v7.3')

