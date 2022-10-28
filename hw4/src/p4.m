clear 
clc
close all

sig = generate_signal(3);
sigL = length(sig);
spc = 10;
chips = sigL/spc;

temp = repelem([1,-1,-1,-1,1,-1,1,1],spc);
ca = [];
for i = 1:6
    ca = [ca temp];
end

fs = 1/sigL;
tInt = 0.02;
intSamps = sigL*tInt;
bs = 1;
be = intSamps;

for i = 1:10

Ecode = circshift(ca(1+(i-1)*intSamps:intSamps*i),spc/2);
Lcode = circshift(ca,-spc/2);

E = sig*Ecode';
L = sig*Lcode';

error = 0.5*(E-L)/(E+L);

end
