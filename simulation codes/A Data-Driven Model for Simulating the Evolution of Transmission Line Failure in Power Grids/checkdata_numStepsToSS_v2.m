% Calculating the steps needed and cumulative number of lines in steady state
% Update
% 04/02: Normalize with the number of source failures
% 04/16 : With 10000 iterations data.
% 12/06/2017 : calculating the steps needed and cumulative number of lines in steady state
% 02/11/2018 : Updated with new data

clc
close all
clear numStepsAndLines
% load state variables
%filename = filename
filename = 'CapData02102018IniF5r0.85e0.45Theta0.2.mat';
%filename = 'Data02102018IniF3r0.85e0.45Theta0.1.mat';
%filename = 'OldCapData02102018IniF3r0.85e0.45Theta0.2.mat';
%filename = 'C:\Users\pankaz\Google Drive\Research\Simulation\PowerGrid\MultiFailure\DataUsed\States2017-MIniF=2DG0.85Alpha0.45Beta0.1.mat'; % *
%filename = 'C:\Users\pankaz\Google Drive\Research\Simulation\PowerGrid\MultiFailure\DataUsed\States2017-MIniF=3DG0.98Alpha0.45Beta0.1.mat'; % *
%filename = 'C:\Users\pankaz\Google Drive\Research\Simulation\PowerGrid\MultiFailure\DataUsed\States2016-MIniF=5DG0.99Alpha0.4Beta0.1.mat'; % *
stateLog = load(filename,'States');
stateLog = struct2cell(stateLog);
stateLogMat = cell2mat(stateLog);
powerGridState = stateLogMat(:,8);
% Load capcity fail sequence
capFailSeq = load(filename,'capalogcell');
capFailSeq = struct2cell(capFailSeq);
initFail = capFailSeq{1,1}{1,1}; %Find initial failed lines
numStepsCasStop = 0; % number of steps needed for stopping
numLineFailedStop = numel(initFail); % number of lnie failed when cascading stopped
counter = 0; % counter
capLimit = numel(capFailSeq{1,1});
for capCellIdx=1:1:capLimit
    if (powerGridState(capCellIdx,1) == -1) % if cascade is stopped
        counter = counter + 1;
        numLineFailedStop = stateLogMat(capCellIdx,1);
        numStepsAndLines(counter,:) = [numStepsCasStop numLineFailedStop];
        numStepsCasStop = 0;
    else % if cascade not stopped
        numStepsCasStop = numStepsCasStop + 1;
    end
end
% numStepsAndLines;
figure(1)
%title('Histogram of number of steps when cascade stops')
avgStep = mean(numStepsAndLines(:,1))
h1 = histogram(numStepsAndLines(:,1),'Normalization','probability');
xlabel('Number of time step when cascade stops');
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',15)
xlim([0 12])
grid on
figure(2)
%title('Histogram of cumulative number of lines failed when cascade stops')
h2 = histogram(numStepsAndLines(:,2));
%meanValue2 = h2.Values * [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29]'
xlabel('Total number of failed lines when cascade stops')
set(gca,'Fontname','Helvetica')
set(gca,'Fontsize',15)
grid on