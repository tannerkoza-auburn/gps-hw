%% Fundamentals of GPS - Homework 2
% Author: Tanner Koza

clear 
clc
close all

%% Problem 1
% Least Squares for a single parameter is the average. 

% Part A & B
a = 5; % actual value

n = 250; % number of measurements
m = 1000; % number of simulations

ahat = zeros(1,m); % preallocation
est_err = zeros(1,m);
P = zeros(1,n);
ahat_ = zeros(1,n);
var = zeros(1,n);

for i = 1:n

    for j = 1:m
    
        y = a + randn(i,1);
    
        ahat(j) = mean(y);
    
        est_err(j) = ahat(j) - a;
    end
    
    P(i) = 1/i;
    ahat_(i) = mean(ahat);
    var(i) = std(est_err)^2;
    
end

figure
plot(P)
hold on
plot(var)
title('Estimate Variance vs. # of Measurements')
xlabel('# of Measurements')
ylabel('Variance')
legend('Variance Function', 'Monte Carlo Results')

clearvars

%% Problem 2

x = [0; 1; 2; 3; 4];
y = [0.181; 2.680; 3.467; 3.101; 3.437];

sigma = 0.4;

num_el = length(x);

% Part A
for i = 1:num_el
    
    H = [1 x(i)]; 
 
    est = (H' * H)^-1 * H' * y(i);

end


% Part B

% Part C

% Part D

clearvars
%% Problem 3

% Part A, B, & C

poshat = [0; 0]; % initial position guess
sigma = 0.5; % range standard deviation
var = sigma^2; % range variance

m = 5000; % number of simulations

poshat_ = zeros(2,m);

for i = 1:m

    while true
    
        x0 = poshat(1);
        y0 = poshat(2);
    
        y = [( 25 + var*randn ) - ( x0^2 + y0^2 );
             ( 65 + var*randn ) - ( (x0 - 10)^2 + y0^2 );
             ( 45 + var*randn ) - ( x0^2 + (y0 - 10)^2 );
             ( 85 + var*randn ) - ( (x0 - 10)^2 + (y0 - 10)^2 )];
    
        G = [2*x0,      2*y0;
             2*(x0-10), 2*y0;
             2*x0,      2*(y0-10)
             2*(x0-10), 2*(y0-10)];
     
        dposhat = (G' * G)^-1 * G' * y;
    
        poshat = poshat + dposhat;
    
        if norm(dposhat) < sigma^2
            break
        end
    
    end

    P = sigma^2.*(G' * G)^-1;

    poshat_(:,i) = poshat;

end

est_err = poshat_ - [3; 4];
var = cov(est_err');
% mean_poshat = mean(poshat_, 2);
% 
% figure
% plot(est_err_(1,:))



% Part D
clearvars

%% Problem 4

opts = detectImportOptions('sv_pos_one_epoch.txt');
opts.DataLines = [3 11];
opts.VariableNames = {'SVs','X','Y','Z','Pseudoranges'};
gps_sv = table2array(readtable('sv_pos_one_epoch.txt',opts));

opts.DataLines = [14 19];
opts.VariableNames = {'SOOPs','X','Y','Z','Pseudoranges'};
soop_sv = table2array(readtable('sv_pos_one_epoch.txt', opts));

% Part A
rho = gps_sv(1:4,5);
svPos = gps_sv(1:4,2:4);
sigma = 0.5;

[pos, clock_bias, P, itr] = gnssPosition(rho, svPos, sigma);

lla = ecef2lla(pos')

figure
geoplot(lla(1), lla(2), 'b*')
geobasemap satellite
hold on

clearvars -except gps_sv soop_sv sigma

% Part B
rho = gps_sv(1:9,5);
svPos = gps_sv(1:9,2:4);

[pos, clock_bias, P, itr] = gnssPosition(rho, svPos, sigma);

lla = ecef2lla(pos')

geoplot(lla(1), lla(2), 'r*')
hold on

clearvars -except gps_sv soop_sv sigma clock_bias

% Part C
rho = gps_sv(1:4,5) - clock_bias;
svPos = gps_sv(1:4,2:4);
sigma = 0.707;

[pos, clock_bias, P, itr] = gnssPosition(rho, svPos, sigma);

lla = ecef2lla(pos')

figure
geoplot(lla(1), lla(2), 'g*')
geobasemap satellite
hold on

% Part D
rho = [gps_sv(1:2,5); soop_sv(1:2,5)];
svPos = [gps_sv(1:2,2:4); soop_sv(1:2,2:4)]

[pos, clock_bias, P, itr] = gnssPosition(rho, svPos, sigma);

lla = ecef2lla(pos')

figure
geoplot(lla(1), lla(2), 'c*')
geobasemap satellite
hold on

clearvars -except gps_sv soop_sv sigma

% Part E
rho = gps_sv(1:4,5) - clock_bias;
svPos = gps_sv(1:4,2:4);
sigma = 0.5;

[pos, clock_bias, P, itr] = gnssPosition(rho, svPos, sigma);

lla = ecef2lla(pos')

figure
geoplot(lla(1), lla(2), 'g*')
geobasemap satellite
hold on

%% Problem 5