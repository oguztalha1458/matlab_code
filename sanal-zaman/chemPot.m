function u = chemPot(Phi,x,dx,px,V0,wd,g)
% Kimyasal potansiyel u = <Phi|H|Phi>/<Phi|Phi>
% H = -1/2 d^2/dx^2 + 1/2 x^2 - V0*exp(-x^2/(2*wd^2)) + g*|Phi|^4
n = abs(Phi).^2;
nrm = sum(n)*dx;
if nrm ~= 0
    % Potansiyel ve etkileşim katkısı (x uzayında)
    V = sum((0.5*x.^2-V0*exp(-x.^2./(2*wd^2))).*n + g*n.^3)*dx;
    % Kinetik katkı (Parseval: momentum uzayında tek FFT ile)
    K = 0.5*sum(px.^2.*abs(fftshift(fft(Phi))).^2)*dx/numel(Phi);
    u = (V+K)/nrm;
else
    u = 0;
end
end
