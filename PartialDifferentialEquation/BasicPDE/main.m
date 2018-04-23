%% 基本的一维偏微分方程求解
clc, clear

m = 0;
x = linspace(0, 1, 20); % x取20个点
t = linspace(0, 2, 20); % 时间t取20个点输出

sol = pdepe(m, @pdefun, @ic, @bc, x, t);
u = sol(:, :, 1);        % 取出答案

% 绘图输出
figure(1)
surf(x, t, u)
title('pde数值解')
xlabel('位置x')
ylabel('时间t')
zlabel('数值解u')

% 与解析解比较
figure(2)
surf(x, t, exp(-t)'*sin(pi*x));
title('解析解')
xlabel('位置x')
ylabel('时间t')
zlabel('数值解u')

% 显示特定点上的解(指定x或t)
figure(3)
M = length(t);   % 显示时间终点上的解
xout = linspace(0, 1, 100);   % 输出点的位置
[uout, dudx] = pdeval(m, x, u(M,:), xout);
plot(xout, uout);
title('末时间时各位置下的解');
xlabel('x')
ylabel('u')
