function [u,Et,E_kin,E_trap,E_strong,E_dip] = chemPotdip(Phi,x,dx,px,g,gdd,lamd,eta)
% Kimyasal potansiyel u ve parçacık başına enerjiler
% H = -1/2 d^2/dx^2 + 1/2*lamd^2*|x|^eta + g*|Phi|^4 - gdd*|Phi|^2
n = abs(Phi).^2;
nrm = sum(n)*dx;
if nrm ~= 0
    % Kinetik enerji (Parseval: momentum uzayında tek FFT ile)
    E_kin = 0.5*sum(px.^2.*abs(fftshift(fft(Phi))).^2)*dx/numel(Phi)/nrm;
    E_trap = sum(0.5*lamd^2*abs(x).^eta.*n)*dx/nrm;
    I4 = sum(n.^3)*dx/nrm;   % <|Phi|^4>
    I2 = sum(n.^2)*dx/nrm;   % <|Phi|^2>
    u = E_kin + E_trap + g*I4 - gdd*I2;
    % Enerjide etkileşim terimleri 1/3 ve 1/2 çarpanıyla girer, u'da girmez
    E_strong = g*I4/3;
    E_dip = -gdd*I2/2;
    Et = E_kin + E_trap + E_strong + E_dip;
else
    [u,Et,E_kin,E_trap,E_strong,E_dip] = deal(0);
end
end
