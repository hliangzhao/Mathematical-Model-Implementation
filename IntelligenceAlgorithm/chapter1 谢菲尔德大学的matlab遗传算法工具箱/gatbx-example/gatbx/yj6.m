clc, clear;
%% 离散二次线性系统最优控制问题
%定义遗传算法参数
GGAP=0.8;              %代沟(Generation gap)
XOVR=1;                %交叉率
NVAR=20;               %变量维数
MUTR=1/NVAR;           %变异率
MAXGEN=2000;           %最大遗传代数(Maximum number of generations)
INSR=0.9;              %插入率
SUBPOP=12;             %子代数目
MIGR=0.2;              %迁移率
MIGGEN=20;             %每20代迁移个体
NIND=20;               %个体数目(Number of individuals)
SEL_F='sus';           %选择函数名
XOV_F='recdis';        %重组函数名
MUT_F='mutbga';        %变异函数名
OBJ_F='objlinq';       %目标函数名
FieldDR=feval(OBJ_F,[],1);                         
%Chrom=crtrp(SUBPOP*NIND,FieldDR);                    %创建初始种群
Chrom=crtrp(SUBPOP*NIND,FieldDR);
gen=0;
trace=zeros(MAXGEN,2);                               %遗传算法性能跟踪
ObjV=feval(OBJ_F,Chrom);                             %计算目标函数值
while gen<MAXGEN                                     %代循环
    trace(gen+1,1)=min(ObjV);
    trace(gen+1,2)=mean(ObjV);
    FitnV=ranking(ObjV,[2,0],SUBPOP);                %分配适应度值(Assign fitness values)
    SelCh=select(SEL_F,Chrom,FitnV,GGAP,SUBPOP);                     %选择
    SelCh=recombin(XOV_F,SelCh,XOVR,SUBPOP);                         %重组
    SelCh=mutate(MUT_F,SelCh,FieldDR,[MUTR],SUBPOP);                 %变异
    ObjVOff=feval(OBJ_F,SelCh);                                      %计算子代目标函数值
    [Chrom, ObjV]=reins(Chrom,SelCh,SUBPOP,[1 INSR],ObjV,ObjVOff);   %替代
    gen=gen+1;
    %在子种群之间迁移个体
    if(rem(gen,MIGGEN)==0)
        [Chrom, ObjV]=migrate(Chrom,SUBPOP,[MIGR, 1, 1],ObjV);
    end
end
[Y,I]=min(ObjV);                      
subplot(211);
plot(Chrom(I,:));
hold on;
plot(Chrom(I,:),'.');grid             %最优控制向量分布图
legend('最优控制向量')
subplot(212);
plot(trace(:,1));hold on;
plot(trace(:,2),'-.');grid
legend('解的变化','种群均值的变化');