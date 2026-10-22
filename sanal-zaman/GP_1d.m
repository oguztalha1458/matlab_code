% TG gazının harmonik tuzaktaki taban durumu (sanal zaman, split-step Fourier)
% mu*Phi = [-1/2 d^2/dx^2 + 1/2 x^2 - V0*exp(-x^2/(2*wd^2)) + g*|Phi|^4]*Phi
clc;clear;
N = 20;
V0 = 0;          % ortadaki Gauss çukurunun derinliği (0: yalnızca harmonik tuzak)
wd = 1.0;        % Gauss çukurunun genişliği
g = N^2*pi^2/2;  % TG etkileşim katsayısı
L = 120;
Ngrid = 2^12;
tau = 0.0001;    % sanal zaman adımı
maxit = 1e6;     % en fazla adım sayısı
dx = L/Ngrid;
x = (-Ngrid/2:Ngrid/2-1)*dx;          % FFT için periyodik ızgara (Ngrid nokta)
px = (-Ngrid/2:Ngrid/2-1)*2*pi/L;     % x ızgarasına karşılık gelen momentumlar
lim = 10^-12;                         % mu'daki bağıl değişim bundan küçükse dur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Başlangıç dalga fonksiyonu
Vpot = 0.5*x.^2;   % harmonik tuzak
initflag = input('Initial? Gaussian:1,ThomasFermi: 2,');
Phi = InitialCond(x,Ngrid,g,Vpot,initflag);
Phi0 = Normalise(Phi,dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bir önceki dalga fonksiyonunu dosyadan okutma
%load N_450_V0_50_wd_1_Phi.dat;
%Phi0= N_450_V0_50_wd_1_Phi(:,2)';
%Phi0 = Normalise(Phi0,dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vext = Vpot - V0*exp(-x.^2./(2*wd^2));   % dış potansiyel: tuzak + Gauss çukuru
expK = exp(-tau*0.5*px.^2);              % kinetik adımın çarpanı (momentum uzayında)
Phi = Phi0;
ustart = chemPot(Phi,x,dx,px,V0,wd,g);
du = 1;
it = 0;
% Taban duruma ulaşana kadar sanal zamanda ilerle
while (abs(du)>lim)
    Phi = exp(-tau*(Vext + g*abs(Phi).^4)).*Phi;             % x uzayında potansiyel adımı
    Phi = real(ifft(ifftshift(expK.*fftshift(fft(Phi)))));   % momentum uzayında kinetik adım
    Phi = Normalise(Phi,dx);                                 % normu tekrar 1 yap
    ustop = chemPot(Phi,x,dx,px,V0,wd,g);                    % yeni kimyasal potansiyel
    du = abs((ustop-ustart)/ustart);
    ustart = ustop;
    it = it + 1;
    if mod(it,10000) == 0
        fprintf('%7d. adım: mu = %.8f, du = %.2e\n', it, ustop, du);
    end
    if it >= maxit || isnan(du)
        warning('Yakınsama sağlanamadı (%d adım). tau değerini küçültün.', it);
        break
    end
end
fprintf('Taban durum: mu = %.8f (%d adım)\n', ustop, it);
Veff = Vext + g*abs(Phi).^4;   % etkin potansiyel
dosya = sprintf('N_%d_V0_%g_wd_%g', N, V0, wd);
wavesonuc = [x' Phi'];
csvwrite([dosya '_Phi.dat'],wavesonuc);
dens = abs(Phi).^2;   % yoğunluk (integrali 1)
denssonuc = [x' dens'];
csvwrite([dosya '_dens.dat'],denssonuc)
