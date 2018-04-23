clc, clear;
%% 收获系统最优控制
%定义遗传算法参数
NVAR=20;               %变量维数
RANGE=[0;200];         %变量范围
GGAP=0.8;              %代沟(Generation gap)
XOVR=1;                %交叉率
MUTR=1/NVAR;           %变异率
MAXGEN=500;            %最大遗传代数(Maximum number of generations)
INSR=0.9;              %插入率
SUBPOP=8;              %子种群数
MIGR=0.2;              %迁移率
MIGGEN=20;             %在子种群与迁移之间20代
NIND=20;               %个体数目(Number of individuals)
SEL_F='sus';           %选择函数名
XOV_F='recdis';        %重组函数名
MUT_F='mutbga';        %变异函数名
OBJ_F='objharv';       %目标函数名
FieldDD=rep(RANGE,[1,NVAR]);                         %译码矩阵
gen=0;
trace=zeros(MAXGEN,2);                               %遗传算法性能跟踪
Chrom=crtrp(SUBPOP*NIND,FieldDD);                    %创建初始种群
ObjV=objharv(Chrom);                                 %计算目标函数值
while gen<MAXGEN                                     %代循环
    FitnV=ranking(ObjV,[2 1],SUBPOP);                %分配适应度值(Assign fitness values)
    SelCh=select(SEL_F,Chrom,FitnV,GGAP,SUBPOP);                     %选择
    SelCh=recombin(XOV_F,SelCh,XOVR,SUBPOP);                         %重组
    SelCh=mutate(MUT_F,SelCh,FieldDD,[MUTR],SUBPOP);                 %变异
    ObjVOff=feval(OBJ_F,SelCh);                                      %计算目标函数值
    [Chrom, ObjV]=reins(Chrom,SelCh,SUBPOP,[1 INSR],ObjV,ObjVOff);   %替代
    gen=gen+1;
    [trace(gen,1),I]=min(ObjV);
    trace(gen,2)=mean(ObjV);
    %在子种群之间迁移个体
    if(rem(gen,MIGGEN)==0)
        [Chrom, ObjV]=migrate(Chrom,SUBPOP,[MIGR, 1, 1],ObjV);
    end
end
[Y,I]=min(ObjV);           %输出最优解及其序号，Y为最优解，I为种群的序号
figure(1);plot(Chrom(I,:));
hold on;grid;
plot(Chrom(I,:),'bo')
figure(2);plot(-trace(:,1));
hold on;
plot(-trace(:,2),'-.');
legend('解的变化','种群均值的变化');
xlabel('迭代次数')