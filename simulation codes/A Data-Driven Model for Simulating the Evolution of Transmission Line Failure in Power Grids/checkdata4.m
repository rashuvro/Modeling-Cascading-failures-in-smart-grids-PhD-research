% Update
% 04/02: Normalize with the number of source failures
% 04/16: With 10000 iterations data.
% 02/11/2018 : Updated with new data and used in the paper

clc
%clear
close all
%M = importdata('States2016-MIniF=1DG0.7Alpha0.1Beta0.1.mat');
%whos -States2017test2-MIniF=1DG0.99Alpha0.4Beta0.1.mat % view contents of mat file

% load state variables
%TrueCaps=[20 80 200 500 800 9900]; % Zhuoyao
TrueCaps=[20 60 120 200 332 9900]; % PD
%filename = 'CapData02102018IniF2r0.85e0.45Theta0.2.mat';
filename = 'CapData02102018IniF5r0.85e0.45Theta0.2.mat';
%filename = 'States2017-MIniF=2DG0.85Alpha0.5Beta0.1.mat';
%filename = 'C:\Users\pankaz\Google Drive\Research\Simulation\PowerGrid\MultiFailure\DataUsed\States2017-MIniF=2DG0.85Alpha0.45Beta0.1.mat'; % *

% simulation parameters
% initFailure = 2; % Initial failure
% r = % Power-grid loading level
% e =  % Capacity estimation error
% theta = %Load-shedding constraint level
%stateLogAll = load(filename)
stateLog = load(filename,'States');
stateLog = struct2cell(stateLog);
stateLogMat = cell2mat(stateLog);
powerGridState = stateLogMat(:,8);
% Load capcity fail sequence
capFailSeq = load(filename,'capalogcell');
capFailSeq = struct2cell(capFailSeq);
casStop = 0;
failedCounter = zeros(numel(TrueCaps)+1,numel(TrueCaps)+1); % failure counters
capLimit = numel(capFailSeq{1,1});
for capCellIdx=1:1:capLimit
    initFail = capFailSeq{1,1}{1,capCellIdx}; %Find initial failed lines
    for j = 1:1:numel(initFail) % find index of failed lines
        oldFailIdx(j) = find(TrueCaps==initFail(j)); % source lines
    end
    if (powerGridState(capCellIdx,1) == -1) % if cascade is stopped
        casStop = casStop +1;
        % counting the stopping
        for k = 1:numel(oldFailIdx)
            failedCounter(oldFailIdx(k),numel(TrueCaps)+1) = ...
                failedCounter(oldFailIdx(k),numel(TrueCaps)+1) + 1;
        end
    else % if cascade not stopped
        nextCellIdx = capCellIdx + 1; % go to next cell
        failedLines = capFailSeq{1,1}{1,nextCellIdx}; %Find next failed lines
        for j = 1:1:numel(failedLines) % find index of failed lines
            nextFailedIdx(j) = find(TrueCaps==failedLines(j));
        end
        % counting the failures
        for k = 1:numel(oldFailIdx)
            for l = 1:numel(nextFailedIdx)
                failedCounter(oldFailIdx(k),nextFailedIdx(l)) = ...
                    failedCounter(oldFailIdx(k),nextFailedIdx(l)) + 1/length(oldFailIdx);
            end
        end
    end
    clear oldFailIdx; % clearing first index
    clear nextFailedIdx; % then clear next index
end
casStop;
failedCounter = failedCounter(1:numel(TrueCaps)+1,1:numel(TrueCaps)+1); % discards 9900 line


% failedCounterAll = failedCounter;
% sumFailMatAll = sum(failedCounterAll,2);
% failProbablilityAll = failedCounterAll./repmat(sumFailMatAll,[1,6])
% row Sum
%
figure
disp('# of failed lines of capacity j at next time step following failure of line of capacity i')
disp('rows = i, columns = j (table 1 in paper)')
failedCountInNeed = [failedCounter(1:numel(TrueCaps)-1,1:numel(TrueCaps)-1)]; % excluding 6-th column
disp(failedCountInNeed)
sumFailMat = sum(failedCountInNeed,2);
failProbablility = failedCountInNeed./repmat(sumFailMat,[1,numel(TrueCaps)-1]);
disp('p_i(j) (rows = i, columns = j)')
disp(failProbablility)
%
filename1 = 'failProb.mat';
save(filename1,'failProbablility');
X = [20 60 120 200 332]; % PD
Y = failProbablility;
plot(X,Y(1,:),'-bx',X,Y(2,:),':k+',X,Y(3,:),'md--',X,Y(4,:),'cp:',X,Y(5,:),'ro-','linewidth',1.5)
legend('{\it{i}} = 20 MW','{\it{i}} = 60 MW','{\it{i}} = 120 MW','{\it{i}} = 200 MW','{\it{i}} = 332 MW')
%legend('i = 20 MW','i = 80 MW','i = 200 MW','i = 500 MW','i = 800 MW')
xlabel('{\it j} (MW)');
ylabel('p_{\it{i}}(\it{j})');
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',18)
grid on
hold on
plot(41,0.46,'o','MarkerSize',30,'linewidth',3)
%title('Failure probability j-th capacity line following i-th capacity line fialure')
%}

%{
    % column Sum
    figure(2)
    sumFailMatCol = sum(failedCounter,1); %column sum
    failProbablilityCol = failedCounter./repmat(sumFailMatCol,[5,1]);
    X1 = X';
    Y1 = failProbablilityCol;
    plot(X1,Y1(:,1),'-bx',X1,Y1(:,2),':k+',X1,Y1(:,3),'gd--',X1,Y1(:,4),'cp:',...
    X1,Y1(:,5),'ro-','linewidth',1.5)
    legend('i = 20 MW','i = 80 MW','i = 200 MW','i = 500 MW','i = 800 MW')
    xlabel('j (MW)');
    ylabel('Probability of line failure, P_j(i)');
    set(gca,'Fontname','Helvetica')
    set(gca,'Fontsize',15)
    grid on
    title('Probability of a single i-th line fialure')
%}

