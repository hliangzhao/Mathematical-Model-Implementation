function [pl, ql, pr, qr] = bc(xl, ul, xr, ur, t)
    % 边界条件函数
    pl = ul;
    ql = 0;
    pr = pi*exp(-t);
    qr = 1;
    
end