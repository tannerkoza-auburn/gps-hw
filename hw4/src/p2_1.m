%% Fundamentals of GPS - Homework 4 - Problem 2

clear
clc
close all

%% Part A

% Sampling Information
fS = 1e6; % Sampling Frequency
tS = 1/fS; % Sampling Period

% Generate Signal
sig = generate_signal(1);
tSig = 0:tS:1-tS;
sig =  sig(1:end-1); %sin(100*2*pi.*tSig);
sigLen = length(sig);

% Time Initialization
intP = 100;
tInt = (fS/sigLen)/intP;
numSamps = tInt/tS;
t = 0:tS:tInt-tS;
tLen = length(t);

% PLL Initialization
omegaNom = 2*pi*150; % Nominal Frequency (Hz)

phi = zeros(intP,1);
phi(1) = pi;
omega = zeros(intP,1);
omega(1) = omegaNom;
eSum = zeros(intP,1);
eSumInt = 0;
eHat = zeros(intP,1);
sGenL = zeros(sigLen,1);

bS = 1;
bE = numSamps;

% Controller Gains
bW = 90;
Ki = 5;
Kp = bW*2*0.707;

for i = 1:intP
    
    % Phase Detector
    s = sig(bS:bE); % Reference Signal
    iSig = sin(omega(i).*t + phi(i)); % Generated Signal
    qSig = cos(omega(i).*t + phi(i));
    I = s .* iSig; 
    Q = s.* qSig;
    eSum(i) = atan2( (sum(Q)/numSamps) , (sum(I)/numSamps) );
    eSumDer = (eSum(i) - eSum(max(i-1,1)) ) /tInt;
    eSumInt = ( eSum(i) - eSum(max(i-1,1)) )/tInt;
%     phi(i+1) = phi(i) + rem(omega(i)*tInt,2*pi);

    % Loop Filter
    eHat(i) = Kp*eSum(i) + Ki*eSumInt;
    omega(i+1) = omega(i) + eHat(i);
    phi(i+1) = phi(i) + omega(i+1)*tInt;
    
    % Log Generated Signal
    sGenL(bS:bE) = iSig;

    % Buffer Update
    bS = bS + numSamps;
    bE = bE + numSamps;

end

plot(omega)

figure
plot(eSum)
hold on
% plot(eHat)

figure
plot(sig)
hold on
plot(sGenL)
axis padded
