
clear all;  close all;

addpath(genpath(pwd))

%%  Dictionary Parameters

Pmax            = 50; %The largest period spanned by the NPDs
Dictionary_pool = {'Ramanujan','NaturalBasis','random' };%Type of the dictionary
Dictionary_type = Dictionary_pool{1};


%%  generate signal

t = 0:1:200;
x1c = 2*cos(2*pi*(1/49)*t);

x1m = cos(2*pi*(1/5)*t);

x1 = x1c + x1m;

x2c = 2*cos(2*pi*(1/49)*t + pi/3);

x2m = cos(2*pi*(1/5)*t+ pi/3);

x2 = x2c + x2m;

x3c = 2*cos(2*pi*(1/49)*t + 2*pi/3);

x3m = cos(2*pi*(1/5)*t+ 2*pi/3);

x3 = x3c + x3m;

data_matrix = [x1' x2' x3'];

% Create figure with 3 subplots
figure
subplot(3,1,1)
plot(t, x1)
title('x1')
subplot(3,1,2)
plot(t, x2, 'r')
title('x2')
subplot(3,1,3)
plot(t, x3,'g')
title('x3')



%% LAPIS

opts.Dictionary_type = Dictionary_type;
opts.Pmax            = Pmax;
opts.lambda_1        = 0.1;
opts.lambda_2        = 0.1;
opts.lambda_3        = 1;
opts.rho             = 1e-3;
opts.visual          = 1;
opts.max_iter        = 100;

[Factor,detected_periods,running_time] = LAPIS(data_matrix,opts);

disp(detected_periods)



