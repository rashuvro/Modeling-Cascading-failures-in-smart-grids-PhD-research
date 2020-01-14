% Input parameter table
clc; clear all; close all;
load tail.mat
% IEEE-118 bus system
NumberOfLines=186;
epsilon=0.05;
% Max number of capacities
Capa=[20 80 200 500 800];
C=length(Capa);
L=60;
N=L*length(Capa)*2;
% Assign capacity weights
w=zeros(1,C);
for i=1:C
    w(i)=C+1-i;
end

deltaT=0.1;
% Weights
Wf=0.5; Wcmax=0.5;

% Choose start and end states INDEX
Fi=[2 3 4 5]; Ci=1;
% Fj=10; Cj=2; Sj=2*C*(Fj-1)+2*Cj;
tspan=0:0.05:10;
DeltaT=tspan(2)-tspan(1);
%%%%%%%%%%%%%%%%%%%
ParameterTable1V3;
%%%%%%%%%%%%%%%%%%%
Fend=zeros(length(DGRatio),length(Alpha),length(Fi));
rate=zeros(length(DGRatio),length(Alpha),length(Fi));
vari=zeros(length(DGRatio),length(Alpha),length(Fi));
idx=0;
Dist=zeros(length(DGRatio),length(Alpha),NumberOfLines);
tic
for m=1:length(DGRatio)
    for n=1:length(Alpha)
        a1=ParaSetting{m,n}(1);
        a2=ParaSetting{m,n}(2);
        a3=ParaSetting{m,n}(3);
        a4=ParaSetting{m,n}(4);
        % If start from 1 failure
        for i=1:NumberOfLines-1
            % j==1 for transient, j==2 for stable states
            % find stability probability parametricly first
            % P_{stop}(F_i)
            if i<=floor(a2*NumberOfLines)
                f1=epsilon + a1*( (a2*NumberOfLines-i)/(a2*NumberOfLines) )^4;
            end
            if i>floor(a2*NumberOfLines) && i<=floor(0.5*NumberOfLines)
                f1=epsilon;
            end
            if i>floor(0.5*NumberOfLines)
                f1=min(1, (epsilon + ( (i-0.5*NumberOfLines)/...
                    (NumberOfLines-0.5*NumberOfLines) )^4) );
            end
            for j=1:2
                % For transiant states j==1
                if j==1
                    % If failure starts from smallest capacity
                    for k=1:C
                        % P_{stop}(C^{\max}_i)
                         f2=max(a4, a3*( (Capa(k)-max(Capa))/max(Capa) )^2 );
                        % Weighted average
                        pStable=Wf*f1 + Wcmax*f2;
                        if pStable>1
                            pStable=1;
                        end
%                         [lambdaS lambdaT]=findSumlambda(pStable,deltaT);
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

        [V,D] = eig(Q);
        U=diag(V);
        dist=zeros(1,NumberOfLines);

        for k=1:length(Fi)
            i=2*C*(Fi(k)-1)+2*Ci-1;
            Fj=Fi(k); Cj=Ci;
            for j=i+1:2:2*NumberOfLines*C
                if Cj==length(Capa)+1
                    Fj=Fj+1;
                    Cj=1;
                end
                if U(j)~=0
                    temp=V(i,j)/U(j);
                    dist(Fj)=dist(Fj)+temp;
                end
                Cj=Cj+1;
            end
        %     dist=dist(dist~=0);
            Dist(m,n,:)=dist;

            for i=Fi(k):NumberOfLines
                Fend(m,n,k)=Fend(m,n,k)+Dist(m,n,i)*i;
            end

            rate(m,n,k)=(Fend(m,n,k)-Fi(k))/Fi(k);
        end
        clear Q
    end
end
toc
save ('RateCommbineV1.mat','rate','DGRatio','Alpha','-v7.3')