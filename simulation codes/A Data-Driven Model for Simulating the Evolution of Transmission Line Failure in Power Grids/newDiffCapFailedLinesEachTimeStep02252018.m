% Ploting average cumulative line failures with different capacities

% Update
% 04/02: Normalize with the number of source failures
% 04/16 : With 10000 iterations data.
% 12/18/2017 : calculating number of failed lines in each catagory at steady state
% 01/11/2018 : calculating number of failed lines of each type at every time steps
% 02/11/2018 : new data with 179 lines
% 02/11/2018 : checking many cases of cascading failure from OPF data

clc;clear;close all
% load state variables
filename = 'CapData02102018IniF4r0.85e0.45Theta0.2.mat';

% simulation parameters
stateLog = load(filename,'States');
stateLog = struct2cell(stateLog);
stateLogMat = cell2mat(stateLog);
powerGridState = stateLogMat(:,8);
% Load capcity fail sequence
capFailSeq = load(filename,'capalogcell');
capFailSeq = struct2cell(capFailSeq);

TrueCaps=[20 60 120 200 332 9900]; % PD
stopPosition = find(powerGridState == -1); % indices of all cascade stop steps
A = stopPosition; % indices of all cascade stop steps
B = diff(A); % time steps of all cascading failurs, ignoring first one
B1 = find(B>0); % taking indices of those cascading process with time step>4

idxtrack = []; % initialize all indices
for ii = 1:length(B1)
    iniFailIdx = A(B1(ii,1),1)+1; % start index of (B1(ii,1)+1)th cascading process
    % taking those cases with initial failures with line capacity greater
    % than 120 and less than 9900
    if sum(capFailSeq{1,1}{1,iniFailIdx} >= 20)>=1 %&& sum(capFailSeq{1,1}{1,iniFailIdx} == 9900)<1
        idxtrack = [idxtrack ii];
        %capFailSeq{1,1}{1,iniFailIdx}
    end
end

% Initialization
s = 1;
stopIdx  = 1;
d = stopPosition(stopIdx,1);

newFailLine1AtTimeSteps = zeros(numel(idxtrack),max(B)); % new failure counters   
newFailLine2AtTimeSteps = zeros(numel(idxtrack),max(B)); % new failure counters
newFailLine3AtTimeSteps = zeros(numel(idxtrack),max(B)); % new failure counters
newFailLine4AtTimeSteps = zeros(numel(idxtrack),max(B)); % new failure counters
newFailLine5AtTimeSteps = zeros(numel(idxtrack),max(B)); % new failure counters
cascadeIdx = 1;
for casStop = idxtrack
    s = A(B1(casStop,1)); % stopping index of the previours cascading failure, 5135 is good choice for iniF=5 data
    d = find(A==s) + 1;   % index of A for the current cascading failure
    timeSteps = A(d,1) - s; % time steps of the cascading failure process, A(d,1) = stoping time step of current cascading failure
    failedLineatTimeSteps = zeros(numel(TrueCaps),timeSteps); % failure counters
    for timeIdx = 1:1:timeSteps
        lineFailed = capFailSeq{1,1}{1,s+timeIdx}; %Find failed lines at time index, s= stoping time index of the previous cascadig failure
        for j  = 1:length(lineFailed)
            c = lineFailed(j,1);
            if c == 20
                failedLineatTimeSteps (1,timeIdx) = failedLineatTimeSteps (1,timeIdx) + 1;
            elseif c == 60
                failedLineatTimeSteps (2,timeIdx) = failedLineatTimeSteps (2,timeIdx) + 1;
            elseif c == 120
                failedLineatTimeSteps (3,timeIdx) = failedLineatTimeSteps (3,timeIdx) + 1;
            elseif c == 200
                failedLineatTimeSteps (4,timeIdx) = failedLineatTimeSteps (4,timeIdx) + 1;
            elseif c == 332
                failedLineatTimeSteps (5,timeIdx) = failedLineatTimeSteps (5,timeIdx) + 1;
            elseif c == 9900
                failedLineatTimeSteps (6,timeIdx) = failedLineatTimeSteps (6,timeIdx) + 1;
            end
        end
    end
      % New line failures of each type at each time step
%     zeroPaddingSize = max(B)-timeSteps; % zero padding
%     newFailLine1AtTimeSteps(cascadeIdx,:)=[failedLineatTimeSteps(1,:) zeros(1,zeroPaddingSize)];
%     newFailLine2AtTimeSteps(cascadeIdx,:)=[failedLineatTimeSteps(2,:) zeros(1,zeroPaddingSize)];
%     newFailLine3AtTimeSteps(cascadeIdx,:)=[failedLineatTimeSteps(3,:) zeros(1,zeroPaddingSize)];
%     newFailLine4AtTimeSteps(cascadeIdx,:)=[failedLineatTimeSteps(4,:) zeros(1,zeroPaddingSize)];
%     newFailLine5AtTimeSteps(cascadeIdx,:)=[failedLineatTimeSteps(5,:) zeros(1,zeroPaddingSize)];

    % Cumulative line failures of each type at each time step
    zeroPaddingSize = max(B)-timeSteps; % zero padding
    cumsum1 = cumsum(failedLineatTimeSteps(1,:));
    cumsum2 = cumsum(failedLineatTimeSteps(2,:));
    cumsum3 = cumsum(failedLineatTimeSteps(3,:));
    cumsum4 = cumsum(failedLineatTimeSteps(4,:));
    cumsum5 = cumsum(failedLineatTimeSteps(5,:));
    cumFailLine1AtTimeSteps(cascadeIdx,:)=[cumsum1 cumsum1(1,end)*ones(1,zeroPaddingSize)];
    cumFailLine2AtTimeSteps(cascadeIdx,:)=[cumsum2 cumsum2(1,end)*ones(1,zeroPaddingSize)];
    cumFailLine3AtTimeSteps(cascadeIdx,:)=[cumsum3 cumsum3(1,end)*ones(1,zeroPaddingSize)];
    cumFailLine4AtTimeSteps(cascadeIdx,:)=[cumsum4 cumsum4(1,end)*ones(1,zeroPaddingSize)];
    cumFailLine5AtTimeSteps(cascadeIdx,:)=[cumsum5 cumsum5(1,end)*ones(1,zeroPaddingSize)];
        
    s  = d;
    stopIdx = stopIdx + 1;
    d = stopPosition(stopIdx,1);
    clear failedLineatTimeSteps; % clearing first index
    cascadeIdx = cascadeIdx+1;
end

% plotting new line failures
%{
figure(1)
subplot(3,2,1)
plot(newFailLine1AtTimeSteps','k-','linewidth',1)
hold on
plot(mean(newFailLine1AtTimeSteps),'r-','linewidth',6)
plot(std(newFailLine1AtTimeSteps),'b:','linewidth',6)
xlabel('Time steps')
ylabel('# of 20 MW line failures')
subplot(3,2,2)
plot(newFailLine2AtTimeSteps','k-','linewidth',1)
hold on
plot(mean(newFailLine2AtTimeSteps),'r-','linewidth',6)
plot(std(newFailLine2AtTimeSteps),'b:','linewidth',6)
xlabel('Time steps')
ylabel('# of 60 MW line failures')
subplot(3,2,3)
plot(newFailLine3AtTimeSteps','k-','linewidth',1)
hold on
plot(mean(newFailLine3AtTimeSteps),'r-','linewidth',6)
plot(std(newFailLine3AtTimeSteps),'b:','linewidth',6)
xlabel('Time steps')
ylabel('# of 120 MW line failures')
subplot(3,2,4)
plot(newFailLine4AtTimeSteps','k-','linewidth',1)
hold on
plot(mean(newFailLine4AtTimeSteps),'r-','linewidth',6)
plot(std(newFailLine4AtTimeSteps),'b:','linewidth',6)
xlabel('Time steps')
ylabel('# of 200 MW line failures')
subplot(3,2,5)
plot(newFailLine5AtTimeSteps','k-','linewidth',1)
hold on
plot(mean(newFailLine5AtTimeSteps),'r-','linewidth',6)
plot(std(newFailLine5AtTimeSteps),'b:','linewidth',6)
xlabel('Time steps')
ylabel('# of 332 MW line failures')
%}

plot(mean(cumFailLine1AtTimeSteps),'-kd','linewidth',2)
hold on
plot(mean(cumFailLine2AtTimeSteps),'--+b','linewidth',2)
plot(mean(cumFailLine3AtTimeSteps),':rs','linewidth',2)
plot(mean(cumFailLine4AtTimeSteps),'-.co','linewidth',2)
plot(mean(cumFailLine5AtTimeSteps),'-m','linewidth',2)
xlabel('Time steps')
ylabel('Total number of failed lines')
legend('20 MW','60 MW','120 MW','200 MW','332 MW')
xlim([1 14])
grid on
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',18)
%plot(time,X(:,1),'-kd',time,X(:,2),'--+b',time,X(:,3),':rs',time,X(:,4),'-.co',time,X(:,5),'-m','linewidth',2)

figure(2)
x = 1:14; 
y1 = mean(cumFailLine1AtTimeSteps);
err1 = std(cumFailLine1AtTimeSteps);
errorbar(x,y1,err1,'-kd','linewidth',2,'MarkerSize',7,...
    'MarkerEdgeColor','black','MarkerFaceColor','black')
hold on
x = 1:14; 
y2 = mean(cumFailLine2AtTimeSteps);
err2 = std(cumFailLine2AtTimeSteps);
errorbar(x,y2,err2,'--+b','linewidth',2,'MarkerSize',7,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue')
x = 1:14; 
y3 = mean(cumFailLine3AtTimeSteps);
err3 = std(cumFailLine3AtTimeSteps);
errorbar(x,y3,err3,':rs','linewidth',2,'MarkerSize',7,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
x = 1:14; 
y4 = mean(cumFailLine4AtTimeSteps);
err4 = std(cumFailLine4AtTimeSteps);
errorbar(x,y4,err4,'-.co','linewidth',2,'MarkerSize',7,...
    'MarkerEdgeColor','cyan','MarkerFaceColor','cyan')

x = 1:14; 
y5 = mean(cumFailLine5AtTimeSteps);
err5 = std(cumFailLine5AtTimeSteps);
errorbar(x,y5,err5,'-m','linewidth',2,'MarkerSize',7,...
    'MarkerEdgeColor','magenta','MarkerFaceColor','magenta')

xlabel('Time steps')
ylabel('Total number of failed lines')
legend('20 MW','60 MW','120 MW','200 MW','332 MW')
xlim([1 14])
grid on
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',18)




