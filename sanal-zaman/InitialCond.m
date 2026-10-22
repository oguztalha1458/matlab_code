function Phi = InitialCond(x,Ngrid,g,Vpot,initflag)
% Başlangıç dalga fonksiyonu: 1 = Gauss, 2 = Thomas-Fermi
if initflag==1
    % Harmonik osilatörün taban durumu (p0112 denk. 3.11, deneme dalga fonksiyonu)
    Phi = (1/pi)^(1/4)*exp(-0.5*x.^2);
elseif initflag==2
    % TG Thomas-Fermi: g*|Phi|^4 = u - Vpot  =>  Phi = ((u - Vpot)/g)^(1/4)
    % Normalizasyon u = N verir; g = N^2*pi^2/2 olduğundan N = sqrt(2*g)/pi.
    u = sqrt(2*g)/pi;                    % kimyasal potansiyel
    Phi = (max(u - Vpot, 0)/g).^(1/4);   % TF yarıçapının dışında sıfır
else
    error('initflag 1 (Gauss) ya da 2 (Thomas-Fermi) olmalı.');
end
end
