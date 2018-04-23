%% 二维状态空间上的偏微分方程求解(亦可使用pdetool工具箱来处理！掌握)
% 正方形区域上的热传导方程
clc, clear

% 问题定义
g = 'squareg';   % 定义正方形区域
b = 'squareb1';  % 边界上为0条件
c = 1; a = 0; f = 0; d = 1;

% 产生初始的三角形网格
[p, e, t] = initmesh(g);

% 定义初始条件
u0 = zeros(size(p, 2), 1);
ix = find(sqrt(p(1,:).^2 + p(2,:).^2) < 0.4);
u0(ix) = 1;

% 在时间段0~0.1上求解
nframe = 20;
tlist = linspace(0, 0.1, nframe);
u1 = parabolic(u0, tlist, b, p, e, t, c, a, f, d);

% 动画显示运算结果
for j = 1:nframe
    pdesurf(p, t, u1(:,j));
    mv(j) = getframe;
end
movie(mv, 10)