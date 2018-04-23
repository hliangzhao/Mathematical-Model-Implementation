%% 使用的solve方法求解微分方程的解析解
clc, clear

syms x y  % 定义符号变量
diff_equ = 'x^2+y+(x-2*y)*Dy=0';
dsolve(diff_equ, 'x')

% 此外，还有常微分方程组，齐次、非齐次常微分方程组
% 用时直接查