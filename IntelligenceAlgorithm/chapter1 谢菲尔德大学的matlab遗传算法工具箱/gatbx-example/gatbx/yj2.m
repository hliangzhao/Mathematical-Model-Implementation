clc, clear
%% 一元单峰函数优化实例
%定义遗传算法参数
NIND=40;               %个体数目(Numbe of individuals)
MAXGEN=500;            %最大遗传代数(Maximum number of generations)
NVAR=20;               %变量的维数
PRECI=20;              %变量的二进制位数(Precision of variables)
GGAP=0.9;              %代沟(Generation gap)
trace=zeros(MAXGEN, 2);
%建立区域描述器(Build field descriptor)
FieldD=[rep([PRECI],[1,NVAR]);rep([-512;512],[1, NVAR]);rep([1;0;1;1],[1,NVAR])];
Chrom=crtbp(NIND, NVAR*PRECI);                       %创建初始种群
gen=0;                                               %代计数器
ObjV=objfun1(bs2rv(Chrom, FieldD));                  %计算初始种群个体的目标函数值
while gen<MAXGEN                                     %迭代
    FitnV=ranking(ObjV);                             %分配适应度值(Assign fitness values)
    SelCh=select('sus', Chrom, FitnV, GGAP);         %选择
    SelCh=recombin('xovsp', SelCh, 0.7);             %重组
    SelCh=mut(SelCh);                                %变异
    ObjVSel=objfun1(bs2rv(SelCh, FieldD));           %计算子代目标函数值 
    [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel);     %重插入
    gen=gen+1;                                                 %代计数器增加
    trace(gen, 1)=min(ObjV);                                   %遗传算法性能跟踪
    trace(gen, 2)=sum(ObjV)/length(ObjV);
end
plot(trace(:,1));hold on;
plot(trace(:,2),'-.');grid;
legend(' 种群均值的变化','解的变化')
%输出最优解及其对应的20个自变量的十进制值,Y为最优解,I为种群的序号
[Y, I]=min(ObjV)
X=bs2rv(Chrom, FieldD);
X(I,:)