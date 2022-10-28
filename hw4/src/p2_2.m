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
% sig = sin(100*2*pi*tSig + pi/3);

% Time Initialization
pInt = 100;
tIntE = sigL/pInt/fS;
nIntSamps = tIntE/tS;
tInt = 0:tS:tIntE-tS;

% PLL Initialization
omega0 = 100*2*pi;
phi0 = deg2rad(0);
bS = 1;
bE = int32(nIntSamps);

% Loop Filter Gains
bW = 1*2*pi;
zeta = 0.9;
Kp = 2*zeta*bW;
Ki = bW^2;

% Preallocation
omega = zeros(pInt,1);
omega(1) = omega0;
phi = zeros(pInt,1);
phi(1) = phi0;
eSum = zeros(pInt,1);
iInt = 0;
qInt = 0;
eSumInt = 0;
eHat = zeros(pInt,1);
iData = zeros(pInt,1);
sGen = zeros(sigL,1);

for i = 1:pInt
    % Phase Detector
    s = sig(bS:bE); % Reference Signal
%     iSig = sin(omega(i).*tInt + phi(i)); % Generated Signal
%     qSig = cos(omega(i).*tInt + phi(i));
    iSig = sin(omega(i).*tSig(bS:bE)); % Generated Signal
    qSig = cos(omega(i).*tSig(bS:bE));

    % Phase Detector Discriminator
    I = s .* iSig; 
    Q = s.* qSig;
    eSum(i) = atan2( (sum(Q)) , (sum(I)) )/(2*pi);
    eSumInt = eSumInt + eSum(i)*tIntE;
   
    % Loop Filter
    eHat(i) = Kp*eSum(i) + Ki*eSumInt;
    omega(i+1) = omega(i) + eHat(i);

    % Extract Data
    iData(i) = sum(iSig);

    % Log Generated Signal
    sGen(bS:bE) = iSig;

    % Phase Bookkeeping
%     phi(i+1) = rem(phi(i) + sin(omega(i)*tIntE),2*pi);

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