function [c, f, s] = pdefun(x, t, u, dudx)
    % 偏微分方程的系数向量函数
    c = pi^2;
    f = dudx;
    s = 0;
end