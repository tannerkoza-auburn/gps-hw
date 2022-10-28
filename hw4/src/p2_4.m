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
bW = 10*2*pi;
zeta = 0.707;
Kp = 2*zeta*bW;
Kd = bW^2;

% Preallocation
omega = zeros(pInt,1);
omega(1) = omega0;
phi = zeros(pInt,1);
phi(1) = phi0;
eSum = zeros(pInt,1);
eSumInt = 0;
eHat = 0;
iData = zeros(pInt,1);
sGen = zeros(sigL,1);
omegaNom = 100;
omega = omegaNom;
phase = 0;
eSumOld = 1.0472;


for i = 1:pInt
    tVec = 0:tS:(tIntE-tS);

    estSin = sin(2*pi*omega*tVec + phase);
    estCos = cos(2*pi*omega*tVec + phase);
    phase = phase + rem(2*pi*omega*tIntE,2*pi);

    phaseVec(i) = phase;

    s = sig(bS:bE);
    I = s .* estSin;
    Q = s.* estCos;

    eSum = atan(sum(Q)/sum(I));

    iVec(i) = sum(I);
    errVec(i) = eSum;

    out = Kp*eSum + Kd*(eSum-eSumOld)/tIntE;
    eSumOld = eSum;
    eHat = eHat + out*tIntE;

    omega = omegaNom + eHat;
    omegaNom = omega;
    omegaVec(i) = omega;

    bS = int32(bS + nIntSamps);
    bE = int32(bE + nIntSamps);
end
