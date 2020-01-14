% calculating number of failed lines in each catagory at steady state and
% Saved data from here will be used for the PYTHON plot

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
casStop = 1;
TrueCaps=[20 60 120 200 332 9900]; % PD
totalStops = sum(powerGridState == -1);
failedCounter = zeros(numel(TrueCaps)+1,totalStops); % failure counters
capLimit = numel(capFailSeq{1,1}); % how many batch/group of capacity failures
for capCellIdx=1:1:capLimit % for each batch/group line failures 
    initFail = capFailSeq{1,1}{1,capCellIdx}; % Find initial failed lines
    for j = 1:1:numel(initFail) % find index of failed lines
        oldFailIdx(j) = find(TrueCaps==initFail(j)); % source lines
    end
    if (powerGridState(capCellIdx,1) == -1) % if cascade is stopped      
        % counting the stopping through casStop variable
        for k = 1:numel(oldFailIdx)
            failedCounter(oldFailIdx(k),casStop) = ...
                failedCounter(oldFailIdx(k),casStop) + 1;
        end
        % Storing total number of failures at each cascade stop at last row
        % of failedCounter matrix
        failedCounter(numel(TrueCaps)+1,casStop) =  sum(failedCounter(:,casStop));       
        casStop = casStop +1; % go to next cascade stop
    else % if cascade not stopped
        failedIdx = oldFailIdx;
        % counting the failures
        for k = 1:numel(failedIdx)
            failedCounter(failedIdx(k),casStop) = ...
                failedCounter(failedIdx(k),casStop) + 1;
        end
        %failedCounter
    end
    clear oldFailIdx; % clearing first index
end

failedCountInNeed = [failedCounter(1:5,:);sum(failedCounter(1:5,:))]';
filename = 'numFailedType.xlsx';
xlswrite(filename,failedCountInNeed)


% find maximum number of lines in all types
%{
max(failedCounter(1,:)) % maximum number of 20 MW failed line among all runs
max(failedCounter(2,:)) % maximum number of 60 MW failed line among all runs
max(failedCounter(3,:)) % maximum number of 120 MW failed line among all runs
max(failedCounter(4,:)) % maximum number of 200 MW failed line among all runs
max(failedCounter(5,:)) % maximum number of 332 MW failed line among all runs
%}
