clc, clear;
%% 多元多峰函数的优化实例:shubert函数
[x1,x2]=meshgrid(-10:.1:10);
figure(1);mesh(x1,x2,shubert(x1,x2));            %画出shubert函数图像
%定义遗传算法参数
NIND=40;               %个体数目(Number of individuals)
MAXGEN=50;             %最大遗传代数(Maximum number of generations)
NVAR=2;                %变量数目
PRECI=25;              %变量的二进制位数(Precision of variables)
GGAP=0.9;              %代沟(Generation gap)
%建立区域描述器(Build field descriptor)
FieldD=[rep([PRECI],[1,NVAR]);rep([-10;10],[1,NVAR]);rep([1;0;1;1],[1,NVAR])];
Chrom=crtbp(NIND, NVAR*PRECI);                         %创建初始种群
gen=0;                                                 
trace=zeros(MAXGEN, 2);                                %遗传算法性能跟踪初始值
x=bs2rv(Chrom, FieldD);                                %初始种群十进制转换
ObjV=shubert(x(:,1),x(:,2));                           %计算初始种群的目标函数值
while gen<MAXGEN
    FitnV=ranking(ObjV);                               %分配适应度值(Assign fitness values)
    SelCh=select('sus',Chrom,FitnV,GGAP);              %选择
    SelCh=recombin('xovsp',SelCh,0.7);                 %重组
    SelCh=mut(SelCh);                                  %变异
    x=bs2rv(SelCh,FieldD);                             %子代十进制转换
    ObjVSel=shubert(x(:,1),x(:,2));
    [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);  %重插入
    gen=gen+1;
    [Y, I]=min(ObjV);
    Y,bs2rv(Chrom(I,:),FieldD)                         %输出每一次的最优解及其对应的自变量值
    trace(gen,1)=min(ObjV);                            %遗传算法性能跟踪
    trace(gen,2)=sum(ObjV)/length(ObjV);
    if(gen==50)                                        %迭代数为50时画出目标函数值分布图
        figure(2);
        plot(ObjV);hold on;
        plot(ObjV,'b*');grid;
    end
end
figure(3);clf;
plot(trace(:,1));hold on;
plot(trace(:,2),'-.');grid
legend('解的变化','种群均值的变化')