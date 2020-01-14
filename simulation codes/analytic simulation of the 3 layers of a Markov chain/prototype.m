clc; clear; close all;

NumberOfLines = 186;
% given r, e and theta
r = 0.95; % r in [0.5, 0.95]
e = 0.25; % e in [0.1. 0.25]
theta = 0.25; % theta in [0, 0.25]

% a's are calculated from r, e and theta
a1 = 0.4 -  0.25 * r - e * (0.2 - e) -  0.25 * theta;
a2 = 0.1 * (1 -  0.5 * r + e * (0.2 - e) -  0.7 * theta);
% epsilon
epsilon = 0.6 - 0.4* r - 0.5 * e - 0.3 * theta;
b = 0.6;

% formulate p(F_i), which is the pStop without human error
p = ones(1, NumberOfLines);
for i = 1 : NumberOfLines-1
    % find the pStop probability f1 first, i.e., with no human error
    if i<=floor(a2*NumberOfLines)
            p(i) = epsilon + a1*( (a2*NumberOfLines-i)/(a2*NumberOfLines) )^2;
    end
    if i>floor(a2*NumberOfLines) && i<=floor(b*NumberOfLines)
            p(i) = epsilon;
    end
    if i>floor(b*NumberOfLines)
            p(i) = min(1, (epsilon + ( (i-b*NumberOfLines)/(NumberOfLines-b*NumberOfLines) )^3) );
    end
end
% figure (1)
% plot(1: NumberOfLines, p)

IniF = 2;

dist = zeros(1, NumberOfLines);
tic
ITER = 100000;
counter = 0;
% pre-allocate States to incredibly accelerate the simulation...
States = zeros(10000000, 4);

for iter = 1 : ITER
    if counter > 0.9999* length(States) % to avoid too many iterations
        break;
    end
    % Specify initial failures
    F = IniF;
    % initial communication factor
    d = 1;
    Y = 3;
    % initial human error probability
    h = 0.01;

    while true
        % determine the pStop for current step and check cascade
        pStop = p(F)*d*(1-h);
        if rand < pStop
            % cascade stop, get out of one iteration/realization
            break;
        else
            % one more transmission-line failure
            F = F + 1;
            % update comm.
            [d, Y] = dfunc(F, d, Y, h);
            % update human
            h = hfunc(F, d, h);
        end
        counter = counter + 1;
        States(counter, 1) = F; % Failures of transmission lines
        States(counter, 2) = Y; % Failures of comm. components
        States(counter, 3) = h; % HEP        
    end
    States(counter, 4) = -1; % cascade stop
    % store the blackout size after one iteration
    dist(F) = dist(F) + 1;
end
toc
dist = dist / ITER;
figure (1)
plot(IniF : NumberOfLines, dist(IniF : NumberOfLines))
figure (2)
loglog(IniF : NumberOfLines, dist(IniF : NumberOfLines))
