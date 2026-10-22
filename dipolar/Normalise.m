function p = Normalise(Phi,dx)
% Dalga fonksiyonunu normalize eder: sum(|Phi|^2)*dx = 1
nrm = sum(abs(Phi).^2)*dx;
if nrm ~= 0
    p = Phi/sqrt(nrm);
else
    p = Phi;
end
end
