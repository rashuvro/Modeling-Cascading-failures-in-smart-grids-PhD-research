% DGRatio=[0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.83 0.86 0.89];
DGRatio=[0.01 0.5 0.7 0.9 0.99];
% Alpha=[0.25 0.35 0.45 0.495];
Alpha=[0.25 0.35 0.45 0.5];
ParaSetting=cell(length(DGRatio),length(Alpha));

% 0.01
ParaSetting{1,1}=[4876 16 4733+4553 3559];
ParaSetting{1,2}=[4669 15.5 4222+4822 3406];
ParaSetting{1,3}=[4054 15 4535+4182 3379];
ParaSetting{1,4}=[3967 14.5 4211+4114 3140];

% 0.5
ParaSetting{2,1}=[4491 14 3907+4065 2650];
ParaSetting{2,2}=[4074 11 4099+3693 2195];
ParaSetting{2,3}=[3186 8 3776+2932 1495];
ParaSetting{2,4}=[235 2 592+510 325];

% 0.7
ParaSetting{3,1}=[4439 12.5 4522+3305 2505];
ParaSetting{3,2}=[3941 10 4065+3520 1952];
ParaSetting{3,3}=[2656 8 3466+2673 1025];
ParaSetting{3,4}=[131 2 368+430 293];

% 0.9
ParaSetting{4,1}=[4215 12 3860+3815 2418];
ParaSetting{4,2}=[3594 10 3721+3710 1764];
ParaSetting{4,3}=[2037 7 1399+1533 715];
ParaSetting{4,4}=[91 2 289+399 288];

% 0.99
ParaSetting{5,1}=[4145 12 3908+3541 2167];
ParaSetting{5,2}=[3125 10 3866+3414 1515];
ParaSetting{5,3}=[1939 7 1333+1441 651];
ParaSetting{5,4}=[81 2 301+383 256];


for i=1:length(DGRatio)
    for j=1:length(Alpha)
        ParaSetting{i,j}(1)=ParaSetting{i,j}(1)/10000;
        ParaSetting{i,j}(2)=ParaSetting{i,j}(2)/186;
        ParaSetting{i,j}(3)=ParaSetting{i,j}(3)/20000;
        ParaSetting{i,j}(4)=ParaSetting{i,j}(4)/10000;
    end
end