%% 刚性微分方程的求解
% 高阶微分方程的解法(难点在于如何转换)
clc, clear

[T, Y] = ode45('example_func2', [0, 1], [0; 1; -1])

% 绘制函数图像
plot(T, Y(:, 1), '-', T, Y(:, 2), '--')
xlabel('time t');
ylabel('solution y');
legend('y1', 'y2');