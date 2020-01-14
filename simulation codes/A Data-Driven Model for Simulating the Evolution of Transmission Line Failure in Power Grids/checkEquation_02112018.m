%pankaz.ece@gmail.com
clc; clear; close all
% using our model
% Time-evolution of failure probabilities of transmission lines 
% Time-evolution of failures of transmission lines at discrete time steps
% Total line failures data from proposed model with different kappas
% 02/11/2018 : new data with 179 lines

% % Loading Probability Vector
filename = 'failProb.mat';
pijAll = load(filename);
pijAll = struct2cell(pijAll);
pijAll = cell2mat(pijAll); % conditional failure prob of lines
%pijAll(5,4) = 0.001;

totTimeStep = 14; % Total time step
N_0 = [38 58 56 20 7]; % num of lines of each capacities (20 60 120 200 332)
N_final = [8.5015 3.597 1.6476 0.3647 0.1265]; % num of lines failed in cascading failure simulation Matpower
iniFailIdx = [1 2 3 4 5]; % index of intitial failure: same as in Matpower OPF
iniFailLine = [0.6551 1.0258 0.9526 0.3447 0.1229]; % number of initial failed line of corresponding index
kappa = 3.55; % time steps when casedade stops
%for kappa = [3.55 8 15 20]; % time steps when casedade stops
%
clear X funLine;
funLine(1,:) = N_final; % Number of functional line at time 0 for all types of CAPs
funLine(1,iniFailIdx)= funLine(1,iniFailIdx) - iniFailLine; % subtracting the failed line (same intitial failure as in Matpower OPF)

% first time instant
fi(1,:) = zeros(1,5); % Failed line is initialized at zero at time 1
fi(1,iniFailIdx) = iniFailLine; % Failed line at time 1, type iniFail
X(1,:) = zeros(1,5); % Fraction of failed lines of each type: Initialize
X(1,iniFailIdx) = iniFailLine; % intial failed line, at time instatnt 1 (same intitial failure as in Matpower OPF)

% for 2nd time instant
for k = 2:1:2 % time steps
    for lineTypeJ = 1:1:5 % different capacity lines
        sumTerm = 0;
        for lineTypeI = 1:1:5
            sumTerm = sumTerm + pijAll(lineTypeI,lineTypeJ)*(X(k-1,lineTypeI));
        end
        fi(k,lineTypeJ) = sumTerm;      
        % Sigmoidal function for probability
        if fi(k,lineTypeJ)<0
            pi(k,lineTypeJ) = 0;
%         elseif fi(k,lineTypeJ)<=1
%             pi(k,lineTypeJ) = fi(k,lineTypeJ);
        elseif fi(k,lineTypeJ)<=kappa
            pi(k,lineTypeJ) = fi(k,lineTypeJ)/kappa;
        else
            pi(k,lineTypeJ) = 1;
        end
        X(k,lineTypeJ) = X(k-1,lineTypeJ) + funLine(k-1,lineTypeJ)*pi(k,lineTypeJ); % failed line of type j
        %X(k,lineTypeJ) = min(X(k,lineTypeJ),N_0(1,lineTypeJ)); % adding satuation
        funLine(k,lineTypeJ) = N_final(1,lineTypeJ) - X(k,lineTypeJ);
    end
end

% 3rd time instant to all
for k = 3:1:totTimeStep % time steps
    for lineTypeJ = 1:1:5 % different capacity lines
        sumTerm = 0;
        for lineTypeI = 1:1:5
            sumTerm = sumTerm + pijAll(lineTypeI,lineTypeJ)*(X(k-1,lineTypeI)-X(k-2,lineTypeI));
        end
        fi(k,lineTypeJ) = sumTerm;      
        % Sigmoidal function for probability
        if fi(k,lineTypeJ)<0
            pi(k,lineTypeJ) = 0;
%         elseif fi(k,lineTypeJ)<=1
%             pi(k,lineTypeJ) = fi(k,lineTypeJ);            
        elseif fi(k,lineTypeJ)<=kappa
            pi(k,lineTypeJ) = fi(k,lineTypeJ)/kappa;
        else
            pi(k,lineTypeJ) = 1;
        end
        X(k,lineTypeJ) = X(k-1,lineTypeJ) + funLine(k-1,lineTypeJ)*pi(k,lineTypeJ); % failed line of type j
        %X(k,lineTypeJ) = min(X(k,lineTypeJ),N_0(1,lineTypeJ)); % adding satuation
        funLine(k,lineTypeJ) = N_final(1,lineTypeJ) - X(k,lineTypeJ);
    end
%     funLine
%     X
end
%
figure('Name','Sigmoidal Function','NumberTitle','off');
time = 1:totTimeStep;
plot(time,pi(:,1),'-kd',time,pi(:,2),'--+b',time,pi(:,3),':rs',time,pi(:,4),'-.co',time,pi(:,5),'-m','linewidth',2)
xlim([1 totTimeStep])
xlabel('Time steps')
ylabel('Probability of failure of a line')
legend('20 MW','60 MW','120 MW','200 MW','332 MW')
grid on
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',18)
%title('Time evolution of probability of failure of a line')
figure('Name','Line failure from model','NumberTitle','off');
%X = ceil(X);
time = 1:totTimeStep;
plot(time,X(:,1),'-kd',time,X(:,2),'--+b',time,X(:,3),':rs',time,X(:,4),'-.co',time,X(:,5),'-m','linewidth',2)
xlim([1 totTimeStep])
xlabel('Time steps')
ylabel('Total number of failed lines')
legend('20 MW','60 MW','120 MW','200 MW','332 MW')
grid on
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',18)
%title('Time evolution of failed line')
%}

% for kappa for loop
%{
figure(3)
time = 1:totTimeStep;
plot(time,sum(X,2),'*-','linewidth',2)
hold on
end
xlim([1 totTimeStep])
xlabel('Time steps')
ylabel('Cumulative line failures')
legend('\kappa = 3.55','\kappa = 10','\kappa = 15','\kappa = 20');
grid on
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',20)
%}
