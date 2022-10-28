clear
clc
close all

sig = generate_signal(3);
sigL = length(sig);
spc = 10;
numChips = sigL/spc;

cGen = repelem([1,-1,-1,-1,1,-1,1,1],spc);
ca = [];


fs = 10/sigL;
ts = 1/fs;
sampInt = 80;   

for i = 1:sampInt/length(cGen)
    ca = [ca cGen];
end

Ki = 2;
Kp = 5.33;
eOld = 0;

shiftSamp = 0;
errorInt = 0;
caShift = ca;

for i = 1:sigL/sampInt
    clf
    % shift code w/ controller output
    caShift = circshift(caShift, shiftSamp);

    % find early and late:
    Ecode = circshift(caShift,spc/2);
    Lcode = circshift(caShift,-spc/2);

    % correlation
    sigCorr = sig(1+(i-1)*sampInt:sampInt*i);
    E = sigCorr*Ecode';
    L = sigCorr*Lcode';

    error(i) = 0.5*(E-L)/(E+L);
    shiftSamp = round(5*error(i));
%  

    %shiftSamp = round(shiftSamp + Kp*error(i) + Ki*errorInt);

    figure(1)
    crosscorr(sig,caShift,sampInt-1)

    figure(2)
    stairs(caShift)
     hold on
    stairs(sigCorr)
    
end

    figure(1)
    crosscorr(sig,ca,sampInt-1)

    figure(2)
    stairs(ca)
     hold on
    stairs(sigCorr)