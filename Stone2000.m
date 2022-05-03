function [P_sp,P_mu]=Stone2000(elv,lat,Pslhl,Fsp)
%This function scales production for latitude and altitude

%converting air pressure to altitude
Ps=1013.25; Ts=288.15; zeta=0.0065; gMR=0.0341;
Pz=Ps.*exp((-gMR./zeta).*(log(Ts)-log(Ts-zeta.*elv)));

% Spallogenic production at index for latitudes ()=0, 10, 20, 30, 40, 50, 60-90;
a = [31.8518 34.3699 40.3153 42.0983 56.7733 69.0720 71.8733];
b = [250.3193 258.4759 308.9894 512.6857 649.1343 832.4566 863.1927];
c = [-0.083393 -0.089807 -0.106248 -0.120551 -0.160859 -0.199252 -0.207069];
d = [7.4260e-5 7.9457e-5 9.4508e-5 1.1752e-4 1.5463e-4 1.9391e-4 2.0127e-4];
e = [-2.2397e-8 -2.3697e-8 -2.8234e-8 -3.8809e-8 -5.0330e-8 -6.3653e-8 -6.6043e-8]; 
M = [0.587 0.600 0.678 0.833 0.933 1.000 1.000];

i=round(lat,-1);
if i<10
    n=1;
elseif i<20
    n=2;
elseif i<30
    n=3;
elseif i<40
    n=4;
elseif i<50
    n=5;
elseif i<60
    n=6;
else
    n=7;
end

%Calculating pressure-Dependent scaling sactors for spallation and muomic production
S_sp=a(n)+b(n).*exp(-Pz./150)+c(n).*Pz+d(n).*Pz.^2+e(n).*Pz.^3;
S_mu=M(n).*exp((1013.25-Pz)./242);

Fm=1-Fsp;

P_sp=Pslhl.*S_sp.*Fsp;
P_mu=Pslhl.*S_mu.*Fm;
end