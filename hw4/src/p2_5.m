%% Fundamentals of GPS - Homework 3 - Problem 2
clear
clc
close all

%% Implementation

% Sampling Information
fs = 1e6;
ts = 1/fs;

% Generate Signal
sig = generate_signal(1);
sigL = length(sig);
tEnd = sigL/fs;
tSig = 0:ts:(tEnd-ts);
% sig = sin(100*2*pi*tSig);

% Time
tInt = 0.01;
nSamps = tInt/ts;
t = 0:ts:(tInt-ts);

% PLL
omega0 = 100;
omega = omega0;
phi0 = 0;
eHat = 0;
phi = phi0;
bs = 1;
be = int32(nSamps);
eSumOld = 0;

% Loop Filter Gains
bW = 10;
zeta = 0.707;
Kp = bW^2;
Kd = 2*zeta*bW;

for i = 1:(tEnd/tInt)
    iSig = sin(2*pi*omega*t + phi);
    qSig = cos(2*pi*omega*t + phi);
    phi = phi + rem(2*pi*omega*tInt,2*pi);
    
    s = sig(bs:be);
    I = s.*iSig;
    Q = s.*qSig;
    eSum = atan(sum(Q)/sum(I));

    data(i) = sum(I);
    eSumVec(i) = eSum;

    out = Kp*eSum + (Kd*(eSum-eSumOld))/tInt;
    eHat = eHat + out*tInt;

    omega = omega0 + eHat;
    omegaVec(i) = omega;

    sGen(bs:be) = iSig;

    bs = int32(bs + nSamps);
    be = int32(be + nSamps);

end

plot(sig)
hold on
plot(sGen)

figure 
plot(eSumVec)

figure
plot(omegaVec)
