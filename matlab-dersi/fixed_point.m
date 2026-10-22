% SABİT NOKTA İTERASYONU (FIXED POINT ITERATION)
clc; clearvars;

g = @(x) sin(x) + 8/4.5;         % 8 - 4.5(x - sin x) = 0 denkleminin x = g(x) hali
x = 2.5;                         % başlangıç değeri
tol = 1e-6; max_iter = 100;

disp(' iter          x            |dx|')
for i = 1:max_iter
    x_eski = x;
    x = g(x);
    fprintf('%4d %15.10f %12.2e\n', i, x, abs(x - x_eski))
    if abs(x - x_eski) < tol, break; end
end
fprintf('\nKök = %.6f (%d iterasyon)\n', x, i)
