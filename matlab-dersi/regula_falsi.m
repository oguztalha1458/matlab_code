% REGULA FALSI (DOĞRU ÇİZEREK KÖK BULMA) YÖNTEMİ
clc; clearvars;

f = @(x) 8 - 4.5*(x - sin(x));   % kökü aranan fonksiyon
a = 1; b = 3;                    % kökü içeren aralık
tol = 1e-6; max_iter = 100;

if f(a)*f(b) > 0, error('Bu aralıkta kök yok.'); end

disp(' iter        a            b            m          f(m)')
for i = 1:max_iter
    m = (a*f(b) - b*f(a))/(f(b) - f(a));   % doğrunun x eksenini kestiği nokta
    fprintf('%4d %12.8f %12.8f %12.8f %11.2e\n', i, a, b, m, f(m))
    if abs(f(m)) < tol, break; end
    if f(a)*f(m) < 0, b = m; else, a = m; end   % kökün olduğu tarafı seç
end
fprintf('\nKök = %.6f (%d iterasyon)\n', m, i)
