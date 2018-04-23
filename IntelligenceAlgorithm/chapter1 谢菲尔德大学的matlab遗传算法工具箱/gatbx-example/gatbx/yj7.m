clc, clear;
%% m个地空导弹火力单元对n批空袭目标进行目标分配
%定义遗传算法参数
NIND=40;                    %个体数目(Number of individuals)
MAXGEN=400;                 %最大遗传代数(Maximum number of generations)
GGAP=0.9;                   %代沟(Generation gap)
trace=zeros(MAXGEN,2);      %遗传算法性能跟踪初始值
BaseV=crtbase(15,8);
Chrom=crtbp(NIND, BaseV)+ones(NIND,15);    %初始种群
gen=0;
ObjV=targetalloc(Chrom);                   %计算初始种群函数值
while gen<MAXGEN
    FitnV=ranking(-ObjV);                  %分配适应度值(Assign fitness values)
    SelCh=select('sus',Chrom,FitnV,GGAP);               %选择
    SelCh=recombin('xovsp',SelCh,0.7);                  %重组
    f=rep([1;8],[1,15]);
    SelCh=mutbga(SelCh, f);SelCh=fix(SelCh);            %变异
    ObjVSel=targetalloc(SelCh);                         %计算子代目标函数值
    [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);   %重插入
    gen=gen+1;
    trace(gen,1)=max(ObjV);                             %遗传算法性能跟踪
    trace(gen,2)=sum(ObjV)/length(ObjV);
end
[Y, I]=max(ObjV);Chrom(I,:),Y                           %最优解及其目标函数值
plot(trace(:,1),'-.');hold on;
plot(trace(:,2));grid
legend('解的变化','种群均值的变化')