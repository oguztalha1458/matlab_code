% NEWTON-RAPHSON YÖNTEMİ
clc; clearvars;

f  = @(x) 8 - 4.5*(x - sin(x));  % kökü aranan fonksiyon
df = @(x) -4.5*(1 - cos(x));     % f'nin türevi
x = 2.5;                         % başlangıç değeri
tol = 1e-6; max_iter = 100;

disp(' iter          x             f(x)')
for i = 1:max_iter
    if df(x) == 0, error('Türev sıfır.'); end
    x = x - f(x)/df(x);   % teğetin x eksenini kestiği nokta
    fprintf('%4d %15.10f %12.2e\n', i, x, f(x))
    if abs(f(x)) < tol, break; end
end
fprintf('\nKök = %.6f (%d iterasyon)\n', x, i)
