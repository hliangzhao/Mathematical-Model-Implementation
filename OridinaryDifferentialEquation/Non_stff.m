%% 求解非刚性的常微分方程的方法
clc, clear
[x, y] = ode45('example_func1', [0, 0.5], 1)  % ode45, ode23, ode113