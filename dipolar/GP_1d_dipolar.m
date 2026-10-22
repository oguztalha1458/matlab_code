% Dipolar TG gazının taban durumu (sanal zaman, split-step Fourier)
% mu*Phi = [-1/2 d^2/dx^2 + 1/2*lamd^2*|x|^eta + g*|Phi|^4 - gdd*|Phi|^2]*Phi
% lamd = 0 iken tuzak kalkar: Łebek tezindeki denklem 3.10 (kendi kendine tutunan damlacık)
clc;clear;
N = 20;
gdd0 = 50;           % dipolar etkileşim şiddeti
gdd = gdd0*N;
lamd = 0.1;          % tuzak şiddeti (0: tuzak yok)
eta = 14;            % tuzak üssü, |x|^eta
g = N^2*pi^2/2;      % TG etkileşim katsayısı
L = 50;
Ngrid = 2^12;
tau = 0.0001;        % sanal zaman adımı (N = 150 ve 200 için 1e-5, 1e-4 ıraksayabilir)
maxit = 1e6;         % en fazla adım sayısı
dx = L/Ngrid;
x = (-Ngrid/2:Ngrid/2-1)*dx;          % FFT için periyodik ızgara (Ngrid nokta)
px = (-Ngrid/2:Ngrid/2-1)*2*pi/L;     % x ızgarasına karşılık gelen momentumlar
lim = 10^-12;                         % mu'daki bağıl değişim bundan küçükse dur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diğer başlangıçlar: Gauss ya da kayıtlı bir sonuç
%Phi0 = (1/pi)^(1/4)*exp(-0.5*x.^2);
%load N=500.dat;
%Phi0= N_500(:,2)';
aa = L/2;
Phi0 = sqrt(1/aa)*cos(pi.*x/L);   % başlangıç: kutunun taban durumu
Phi0 = Normalise(Phi0,dx);
Vtrap = 0.5*lamd^2*abs(x).^eta;   % tuzak potansiyeli
expK = exp(-tau*0.5*px.^2);       % kinetik adımın çarpanı (momentum uzayında)
Phi = Phi0;
ustart = chemPotdip(Phi,x,dx,px,g,gdd,lamd,eta);
du = 1;
it = 0;
% Taban duruma ulaşana kadar sanal zamanda ilerle
while (abs(du)>lim)
    n = abs(Phi).^2;
    Phi = exp(-tau*(Vtrap + g*n.^2 - gdd*n)).*Phi;           % x uzayında potansiyel adımı
    Phi = real(ifft(ifftshift(expK.*fftshift(fft(Phi)))));   % momentum uzayında kinetik adım
    Phi = Normalise(Phi,dx);                                 % normu tekrar 1 yap
    % Yeni dalga fonksiyonunun kimyasal potansiyeli ve enerjileri
    [ustop,Et,E_kin,E_trap,E_strong,E_dip] = chemPotdip(Phi,x,dx,px,g,gdd,lamd,eta);
    du = abs((ustop-ustart)/ustart);
    ustart = ustop;
    it = it + 1;
    if mod(it,10000) == 0
        fprintf('%7d. adım: mu = %.6f, du = %.2e\n', it, ustop, du);
    end
    if it >= maxit || isnan(du)
        warning('Yakınsama sağlanamadı (%d adım). tau değerini küçültün.', it);
        break
    end
end
x_rms = sqrt(trapz(x, x.^2.*abs(Phi).^2));   % bulutun genişliği
fprintf('mu = %.4f   E_toplam = %.4f   x_rms = %.4f   (%d adım)\n', ustop, Et, x_rms, it);
fprintf('E_kin = %.4f   E_trap = %.4g   E_strong = %.4f   E_dip = %.4f\n', E_kin, E_trap, E_strong, E_dip);
plot(x,Phi)
xlabel('x'); ylabel('\Phi(x)');
y = [x; Phi];
dosya = sprintf('GP_N_%d_lamda_%s_gdd_%g.dat', N, strrep(num2str(lamd),'.',''), gdd0);
fileID = fopen(dosya,'w');
fprintf(fileID,'%6.4f %10.6f\n',y);
fclose(fileID);
