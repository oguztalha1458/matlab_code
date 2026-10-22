% SECANT (KİRİŞ) YÖNTEMİ
clc; clearvars;

f = @(x) 8 - 4.5*(x - sin(x));   % kökü aranan fonksiyon
x0 = 2.5; x1 = 4;                % iki başlangıç değeri
tol = 1e-6; max_iter = 100;

disp(' iter          x             f(x)')
for i = 1:max_iter
    x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0));   % kirişin x eksenini kestiği nokta
    fprintf('%4d %15.10f %12.2e\n', i, x2, f(x2))
    if abs(f(x2)) < tol, break; end
    x0 = x1; x1 = x2;   % son iki noktayla devam et
end
fprintf('\nKök = %.6f (%d iterasyon)\n', x2, i)
