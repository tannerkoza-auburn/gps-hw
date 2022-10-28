%% Fundamentals of GPS - Homework 3 - Problem 2
clear
clc
close all

%% Part A

% Sampling Information
fS = 1e6;
tS = 1/fS;

% Generate Signal
sig = generate_signal(1);
sigL = length(sig);
tEnd = sigL/fS;
tSig = 0:tS:tEnd-tS;
% sig = sin(102*2*pi*tSig + pi/3);

% Time Initialization
pInt = 200;
tIntE = sigL/pInt/fS;
nIntSamps = tIntE/tS;
tInt = 0:tS:tIntE-tS;

% PLL Initialization
omega0 = 100*2*pi;
phi0 = deg2rad(60);
bS = 1;
bE = int32(nIntSamps);

% Loop Filter Gains
bW = 1*2*pi;
zeta = 0.707;
%Kp = 2*zeta*bW;
%Ki = bW^2;
Kp = 1;
Ki = 5;

% Preallocation
omega = zeros(pInt,1);
omega(1) = omega0;
phi = zeros(pInt,1);
phi(1) = phi0;
eSum = zeros(pInt,1);
eSumInt = 0;
eHat = zeros(pInt,1);
iData = zeros(pInt,1);
sGen = zeros(sigL,1);

for i = 1:pInt
    s = sig(bS:bE);
    iSig = sin(omega(i).*tSig(bS:bE));
    qSig = cos(omega(i).*tSig(bS:bE));

    I = s .* iSig; 
    Q = s.* qSig;
    eSum(i) = atan2(sum(Q),sum(I))/(1);
    eSumInt = eSumInt + eSum(i)*tIntE;
%     eSumDer = (eSum(i) - eSum(max(i-1,1)))/tIntE;

    eHat(i) = 2*(Kp*eSum(i) + Ki*eSumInt);
    omega(i+1) = omega(i) + eHat(i);

    % Log Generated Signal
    sGen(bS:bE) = iSig;

    % Phase Bookkeeping
%     phi(i+1) = phi(i) + rem(omega(i)*tSig(bE),2*pi);

    % Buffer Update
    bS = int32(bS + nIntSamps);
    bE = int32(bE + nIntSamps);

end

figure
plot(eSum)
figure
plot(omega/2/pi)
figure
plot(sig)
hold on
plot(sGen)
axis padded
