% calculating ratio of number of eventually failed lines of each type
% during cascading failures

% Update
% 04/02: Normalize with the number of source failures
% 04/16 : With 10000 iterations data.
% 12/18/2017 : calculating number of failed lines in each catagory at steady state
% 02/11/2018 : Updated with new data and saving for PYTHON plot

clc
close all
% load state variables
filename = 'CapData02102018IniF5r0.85e0.45Theta0.2.mat';
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
cascade = 1; % enumeration of cascade run
newCascade = 0; % tracking new cascade
TrueCaps=[20 60 120 200 332 9900]; % PD
totalStops = sum(powerGridState == -1);
failedCounter = zeros(numel(TrueCaps)+1,totalStops); % failure counters
capLimit = numel(capFailSeq{1,1}); % how many batch/group of capacity failures
initFailedLine = capFailSeq{1,1}{1,1}; % Find initial failed lines
for capCellIdx = 2:1:capLimit % for each batch/group line failures, staring at 2 for skipping initial failure
    if newCascade == 1
        initFailedLine = capFailSeq{1,1}{1,capCellIdx}; % Find initial failed lines
        newCascade = 0;
        % do nothing for this loop, skip the initial capacities
    else
        capCellIdxNew = capCellIdx; 
        if (powerGridState(capCellIdxNew,1) == -1) % if cascade is stopped
            failedLine = capFailSeq{1,1}{1,capCellIdxNew}; % Find initial failed lines
            for j = 1:1:numel(failedLine) % find index of failed lines
                oldFailIdx(j) = find(TrueCaps==failedLine(j)); % source lines
            end
            % counting the stopping for each simulation run through casStop variable
            for k = 1:numel(oldFailIdx)
                failedCounter(oldFailIdx(k),cascade) = ...
                    failedCounter(oldFailIdx(k),cascade) + 1;
            end
            % Storing total number of failures at each cascade stop at last row
            % of failedCounter matrix
            failedCounter(numel(TrueCaps)+1,cascade) =  sum(failedCounter(:,cascade));
            cascade = cascade +1; % go to next cascade run
            newCascade = 1; % for tracking new cascade
        else % if cascade not stopped
            failedLine = capFailSeq{1,1}{1,capCellIdxNew}; % Find initial failed lines
            for j = 1:1:numel(failedLine) % find index of failed lines
                oldFailIdx(j) = find(TrueCaps==failedLine(j)); % source lines
            end
            failedIdx = oldFailIdx;
            % counting the failures
            for k = 1:numel(failedIdx)
                failedCounter(failedIdx(k),cascade) = ...
                    failedCounter(failedIdx(k),cascade) + 1;
            end
            newCascade = 0; % tracking new cascade
            %failedCounter
        end
        clear oldFailIdx; % clearing first index
    end
end
failedCountInNeed = [failedCounter(1:5,:);sum(failedCounter(1:5,:))]';
% find maximum number of lines in all types
max(failedCounter(1,:)); % maximum number of 20 MW failed line among all runs
max(failedCounter(2,:)); % maximum number of 60 MW failed line among all runs
max(failedCounter(3,:)); % maximum number of 120 MW failed line among all runs
max(failedCounter(4,:)); % maximum number of 200 MW failed line among all runs
max(failedCounter(5,:)); % maximum number of 332 MW failed line among all runs

% Ratio
max(failedCounter(1,:))/38 % maximum number of 20 MW failed line among all runs
max(failedCounter(2,:))/58 % maximum number of 60 MW failed line among all runs
max(failedCounter(3,:))/56 % maximum number of 120 MW failed line among all runs
max(failedCounter(4,:))/20 % maximum number of 200 MW failed line among all runs
max(failedCounter(5,:))/7 % maximum number of 332 MW failed line among all runs
