close all; clear; clc;

% solution of y'(t) = 1/(1+y(t)^2) with y(0) = x

syms x t
F(x, t) = (sqrt((81*(x + (x^3)/3) + 81*t)^2 + 2916) ...
    + 81*(x + (x^3)/3) + 81*t)^(1/3)/(3*2^(1/3)) ...
    - (3*2^(1/3))/(sqrt((81*(x + (x^3)/3) + 81*t)^2 + 2916) ...
    + 81*(x + (x^3)/3) + 81*t)^(1/3);
DF = diff(F, x);
DDF = hessian(F, x);

DDF_m = matlabFunction(DDF);
DF_m = matlabFunction(DF);
F_m = matlabFunction(F);

figure();
fsurf(F_m, [-5.0, 5.0, 0.0, 2.0]);
figure();
fsurf(DDF_m, [-5.0, 5.0, 0.0, 2.0]);
xlabel('x');

tmax = 2.0;

[xopt1, fopt1] = fminbnd(@(x) -abs(DDF_m(x, tmax)), -5.0, 5.0)

F1 = @(x, t) 2*DF_m(x, t).^2 + 2.0*DDF_m(x, t).*F_m(x, t);
F2 = @(x, t) integral(@(s) F1(x, s), 0.0, t);
[xopt2, fopt2] = fminbnd(@(x) -abs(F2(x, tmax)), -5.0, 5.0)

F2t = @(x, t) arrayfun(@(t) 2.0*F2(x, t), t);
F3 = @(x, t) integral(@(s) F2t(x, s), 0.0, t);
[xopt3, fopt3] = fminbnd(@(x) -abs(F3(x, tmax)), -5.0, 5.0)