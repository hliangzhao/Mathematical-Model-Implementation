clc, clear;
%% 双积分的优化问题
%定义遗传算法参数
Dim=20;               %变量维数
NIND=20;              %个体数目(Number of individuals)
Preci=20;             %变量的二进制位数(Precision of variables)
MAXGEN=100;            %最大遗传代数(Maximum number of generations)
GGAP=0.8;             %代沟(Generation gap)
SEL_F='sus';          %选择函数名
XOV_F='xovsp';        %重组函数名
MUT_F='mut';          %变异函数名
OBJ_F='objdopi';      %目标函数名
FieldDR=feval(OBJ_F,[],1);                    %计算目标函数值
%建立区域描述器(Build field descriptor)
FieldDD=[rep([Preci],[1,Dim]);FieldDR;rep([1;0;1;1],[1,Dim])];
Chrom=crtbp(NIND, Dim*Preci);                 %创建初始种群
gen=0;
Best=NaN*ones(MAXGEN,1);                      %最优解初值
while gen<MAXGEN                              %最大循环次数
    ObjV=feval(OBJ_F,bs2rv(Chrom,FieldDD));   %计算目标函数值
    Best(gen+1)=min(ObjV);                    %最优解
    plot(log10(Best),'bo');
    FitnV=ranking(ObjV);                      %分配适应度值(Assign fitness values)
    SelCh=select(SEL_F,Chrom,FitnV,GGAP);     %选择
    SelCh=recombin(XOV_F,SelCh);              %重组
    SelCh=mutate(MUT_F,SelCh);                %变异
    Chrom=reins(Chrom,SelCh);                 %重插入
    gen=gen+1;
end
grid;
xlabel('迭代次数');ylabel('目标函数值(取对数)');