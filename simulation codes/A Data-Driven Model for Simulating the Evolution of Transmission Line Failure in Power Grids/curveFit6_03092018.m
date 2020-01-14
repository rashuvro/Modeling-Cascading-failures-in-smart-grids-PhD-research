% plot with power law function
% 02/11/2018 : Updated with new data and used in the paper
clc
clear
close all
filename = 'failProb.mat';
p_ij = load(filename);
stateLog = struct2cell(p_ij);
p_ij = cell2mat(stateLog);
%capacity = [20 80 200 500 800];
capacity = [20 60 120 200 332];
c_j = capacity;
for j_th = 41:41; % Threshold
    c0 = [140 0.5 2.5]';
    for i=1:1:5
        c_i = capacity(1,i);
        fun = @(c) c(1)*c_j.^c(2).*(1 + c(3)*c_i*(c_j - j_th)) - p_ij(i,:);
        %fun = @(c) c(1)*exp(-1*c(2)*c_j) - p_ij(i,:);
        c=lsqnonlin(fun,c0);
        %     plot(c_j,p_ij(i,:),'ro',c_j,c(1)*c_j.^c(2).*(1 + c(3)*c_i*(j_th-c_j)),'b-','linewidth',2)
        %     xlabel('j')
        %     ylabel('fitted function')
        %     hold on
        %     pause(1)
        cavg(i,:)=c';
    end
    %cavg
    %c(1) = 12;
    %c(3) = c(3) + 1.7e-06;
    %
    % grid on
    % legend('Data','Best fit')
    %c = [1.11 0.018 0.08e-6];
    
    %
    figure
    X = [20 60 120 200 332];
    Y = p_ij;
    plot(X,Y(1,:),'bx',X,Y(2,:),'k+',X,Y(3,:),'md',X,Y(4,:),'cp',X,Y(5,:),'ro','linewidth',1.5)
    legend('{\it{i}} = 20 MW','{\it{i}} = 60 MW','{\it{i}} = 120 MW','{\it{i}} = 200 MW','{\it{i}} = 332 MW')
    hold on
    plot(c_j,c(1)*c_j.^c(2).*(1 + c(3)*c_i*(c_j - j_th)),'b-','linewidth',2.5)
    xlabel('{\it j} (MW)');
    ylabel('p_{\it{i}}(\it{j})');
    ylim([0 1])
    l=legend('{\it{i}} = 20 MW','{\it{i}} = 60 MW','{\it{i}} = 120 MW','{\it{i}} = 200 MW','{\it{i}} = 332 MW');
    %title('Power law fitting of line failure probability')
    grid on
    set(gca,'Fontname','Helvetica')
    set(gca,'Fontsize',18)
    %}
    p_ij_estmted = c(1)*c_j.^c(2).*(1 + c(3)*c_i*(j_th-c_j));
    
    MSE = 0;
    for j=1:5
        MSE = MSE + ((norm(p_ij(:,j)-ones(5,1)*p_ij_estmted(1,j)))^2)/5;
    end
    MSE
end
alpha = c(1)
beta = c(2)
gamma = c(3)




