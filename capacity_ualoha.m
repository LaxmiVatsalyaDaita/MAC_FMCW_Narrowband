function Cualoha = capacity_ualoha(Tg,h,fH,lambda)
%Typical values:Tg = 45musec, h = 12MHz/musec, fH = 10MHz

c = 3*10^2;  %m/musec
fc = 77*10^3;  %MHz
Pt = 12; %dBm
Gt = 19.25; %dBi
Gr = 19.25;  %dBi
xit = 20; %m^2

Bc = Tg*h;  %MHz
dtmax = fH*c/(2*h); %m
Deltad = c/(2*Bc);   %m
N = dtmax/Deltad;

dn = ((1:N)-0.5)*Deltad; %m
sigman2 = Pt+Gt+Gr+10*log10(c^2*xit./(64*pi^3*fc^2*dn.^4)); %dBm

drw = 3;  %m
dimax = 1000;  %m  
thetabw = pi/4; %radian

sigma02 = Pt+Gt+Gr+10*log10(c^2*(pi/2-atan(drw/sqrt(dimax^2-drw^2)))/(16*pi^2*fc^2*drw*sqrt(dimax^2-drw^2)));  %dBm

n0 = find(sigman2 == min(sigman2(sigman2 >= sigma02)))
pthtd = 0.999;
%eta = sigman2 + 10*log10(-log(1-pthtd)); %eta used in the paper
eta = sigman2 + 10*((1:N <= n0).*log10(-log(pthtd))+(1:N > n0 & 1:N <= N).*log10(-log(1-pthtd)));  %correct eta

plot(1:N,sigman2,1:N,sigma02*ones(1,N),1:N,eta,'--')

pfan = zeros(1,N);
pspa = (thetabw/(pi-thetabw))^2;
pfan(1:n0) = sqrt(max(0,(10.^((Pt+Gt+Gr-eta(1:n0))/10)*c^2/(16*pi^2*fc^2)-drw^2)/(dimax^2-drw^2)));
pfan(n0+1:N) = 1 - sqrt(min(1,(10.^((Pt+Gt+Gr-eta(n0+1:N))/10)*c^2/(16*pi^2*fc^2)-drw^2)/(dimax^2-drw^2)));
asys = dimax*pspa*(fH/(N*Bc))*sum(pfan);
Cualoha = lambda*exp(-lambda*asys);
end