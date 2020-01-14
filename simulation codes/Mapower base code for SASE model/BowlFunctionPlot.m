clc; clear all; close all;

% load States2016-1OPFhuman600Alpha0.4Beta0.4V1.mat
load States2012-10OPFother990Alpha0.2Beta0.1V1.mat

% For 118 case
NumberOfLines = 186;

total_states = zeros(1, NumberOfLines);
stable_states = zeros(1, NumberOfLines);

for i = 1:length(States)
    total_states(States(i,1)) = total_states(States(i,1)) + 1;
    if States(i,8) == -1
        stable_states(States(i,1)) = stable_states(States(i,1)) + 1;
    end
end

cascade_stop = zeros(1, NumberOfLines);
for i = 1:NumberOfLines
    if total_states(i) ~= 0
        cascade_stop(i) = stable_states(i)/total_states(i);
    end
end

figure (1)
plot(5:NumberOfLines, cascade_stop(5:NumberOfLines))
xlabel('Number of failed transmission lines')
ylabel('Cascade-stop probability')

% clear;
% load States2016-1OPFother600Alpha0.4Beta0.4V3.mat
% 
% % For 118 case
% NumberOfLines = 231;
% 
% total_states = zeros(1, NumberOfLines);
% stable_states = zeros(1, NumberOfLines);
% 
% for i = 1:length(States)
%     total_states(States(i,1)) = total_states(States(i,1)) + 1;
%     if States(i,8) == -1
%         stable_states(States(i,1)) = stable_states(States(i,1)) + 1;
%     end
% end
% 
% cascade_stop = zeros(1, NumberOfLines);
% for i = 1:231
%     if total_states(i) ~= 0
%         cascade_stop(i) = stable_states(i)/total_states(i);
%     end
% end
% 
% figure (2)
% plot(5:230, cascade_stop(5:230))
% xlabel('Number of failed transmission lines')
% ylabel('Cascade-stop probability')
