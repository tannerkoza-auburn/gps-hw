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

% Time
tInt = 0.01;
pInt = floor(tEnd/tInt);
nSamps = tInt/ts;
t = 0:ts:(tInt-ts);

% PLL
omega0 = 99; % Hz
omega = omega0;
phi0 = 0;
phi = phi0;
eSumInt = 0;
eHat = 0;
bs = 1;
be = int32(nSamps);

% Loop Filter Gains
bW = 1*2*pi;
zeta = 0.9;
Kp = 2*zeta*bW;
Ki = bW^2;

% Log Prellocation
phiL = zeros(pInt,1);
eSumL = zeros(pInt,1);
omegaL = zeros(pInt,1);

for i = 1:pInt

    % Oscillator
    iSig = sin(2*pi*omega*t + phi);
    qSig = cos(2*pi*omega*t + phi);
    phi = phi + rem(2*pi*omega*tInt,2*pi);

    phiL(i) = phi; % Phase Log
    
    % Phase Detector Discriminator
    s = sig(bs:be);
    I = s.*iSig;
    Q = s.*qSig;
    eSum = atan(sum(Q)/sum(I));
    eSumInt = eSumInt + eSum*tInt;

    eSumL(i) = eSum; % Phase Error Log

    % Loop Filter
    eHat = Kp*eSum + Ki*eSumInt;
    omega = omega0 + eHat;

    omegaL(i) = omega; % Frequency Log
    
    % Index Update
    bs = int32(bs + nSamps);
    be = int32(be + nSamps);

end

%% Plotting

figure
plot(rad2deg(phiL))
title('Phase')
xlabel('Integration Periods')
ylabel('Phase (degs)')

figure
plot(eSumL)
title('Phase Error')
xlabel('Integration Periods')
ylabel('Phase Error (degs)')

figure
plot(omegaL)
title('In-Phase Frequency')
xlabel('Integration Periods')
ylabel('Frequency (Hz)')

figure
plot(s)
hold on
plot(iSig)
title('Actual vs. In-Phase Signal: Final Integration Period')
xlabel('Time Steps')
ylabel('Amplitude')
legend('Actual','In-Phase')
