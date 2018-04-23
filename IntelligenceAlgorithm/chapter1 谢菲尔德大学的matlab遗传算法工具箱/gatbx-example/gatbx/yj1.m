clc, clear;
%% 简单一元函数优化实例
figure(1);
fplot('variable.*sin(10*pi*variable)+2.0',[-1,2]);   %画出函数曲线
%定义遗传算法参数
NIND=40;        %个体数目(Number of individuals)
MAXGEN=25;      %最大遗传代数(Maximum number of generations)
PRECI=20;       %变量的二进制位数(Precision of variables)
GGAP=0.9;       %代沟(Generation gap)
trace=zeros(2, MAXGEN);                        %寻优结果的初始值
FieldD=[20;-1;2;1;0;1;1];                      %区域描述器(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %初始种群
gen=0;                                         %代计数器
variable=bs2rv(Chrom, FieldD);                 %计算初始种群的十进制转换
ObjV=variable.*sin(10*pi*variable)+2.0;        %计算目标函数值

while gen<MAXGEN
   FitnV=ranking(-ObjV);                                  %分配适应度值(Assign fitness values)         
   SelCh=select('sus', Chrom, FitnV, GGAP);               %选择
   SelCh=recombin('xovsp', SelCh, 0.7);                   %重组
   SelCh=mut(SelCh);                                      %变异
   variable=bs2rv(SelCh, FieldD);                         %子代个体的十进制转换
   ObjVSel=variable.*sin(10*pi*variable)+2.0;             %计算子代的目标函数值
   [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel); %重插入子代的新种群
   variable=bs2rv(Chrom, FieldD);
   gen=gen+1;                                             %代计数器增加
   %输出最优解及其序号，并在目标函数图像中标出，Y为最优解,I为种群的序号
   [Y, I]=max(ObjV);hold on;
   plot(variable(I), Y, 'bo');
   trace(1, gen)=max(ObjV);                               %遗传算法性能跟踪
   trace(2, gen)=sum(ObjV)/length(ObjV);
end
variable=bs2rv(Chrom, FieldD);                            %最优个体的十进制转换
hold on, grid;
plot(variable,ObjV,'b*');
figure(2);
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid
legend('解的变化','种群均值的变化')