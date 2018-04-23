%% 调用MATLAB最优化工具箱函数fmincon进行非线性寻优
function [ ret ] = nonlinear( chrom, sizepop )
for i = 1: sizepop
    x = fmincon(inline('-5*sin(x(1))*sin(x(2))*sin(x(3))*sin(x(4))*sin(x(5))-sin(5*x(1))*sin(5*x(2))*sin(5*x(3))*sin(5*x(4))*sin(5*x(5))+8'), ...
    chrom(i, :)', [], [], [], [], [0 0 0 0 0], [2.8274 2.8274 2.8274 2.8274 2.8274]);
    ret(i, :) = x';
end


end

