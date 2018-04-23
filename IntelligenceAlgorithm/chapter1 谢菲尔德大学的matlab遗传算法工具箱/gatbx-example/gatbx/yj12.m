clc, clear;
%% 多目标优化问题，含有两个优化目标的多目标优化问题
NIND=100;                              %个体数目(Number of individuals)
MAXGEN=50;                             %最大遗传代数(Maximum number of generations)
NVAR=2;                                %变量个数
PRECI=20;                              %变量的二进制位数(Precision of variables)
GGAP=0.9;                              %代沟(Generation gap)
trace1=[];trace2=[];trace3=[];         %性能跟踪
%建立区域描述器(Build field descriptor)
FieldD=[rep([PRECI],[1,NVAR]);[1,1;4,2];rep([1;0;1;1],[1,NVAR])];
Chrom=crtbp(NIND,NVAR*PRECI);          %初始种群
v=bs2rv(Chrom,FieldD);                 %初始种群十进制转换
gen=1;
while gen<MAXGEN
    [NIND, N]=size(Chrom);
    M=fix(NIND/2);
    ObjV1=f1(v(1:M,:));                %分组后第一目标函数值
    FitnV1=ranking(ObjV1);             %分配适应度值(Assign fitness values)
    SelCh1=select('sus',Chrom(1:M,:),FitnV1,GGAP);                 %选择
    ObjV2=f2(v(M+1:NIND,:));           %分组后第二目标函数值
    FitnV2=ranking(ObjV2);
    SelCh2=select('sus',Chrom((M+1):NIND,:),FitnV2,GGAP);          %选择
    SelCh=[SelCh1;SelCh2];             %合并
    SelCh=recombin('xovsp',SelCh,0.7); %重组
    Chrom=mut(SelCh);                  %变异
    v=bs2rv(SelCh,FieldD);
    
    trace1(gen,1)=min(f1(v));
    trace1(gen,2)=sum(f1(v))/length(f1(v));
    trace2(gen,1)=min(f2(v));
    trace2(gen,2)=sum(f2(v))/length(f2(v));
    trace3(gen,1)=min(f1(v)+f2(v));
    trace3(gen,2)=sum(f1(v))/length(f1(v))+sum(f2(v))/length(f2(v));
    gen=gen+1;
end
figure(1);clf;
plot(trace1(:,1));hold on;plot(trace1(:,2),'-.');
plot(trace1(:,1),'.');plot(trace1(:,2),'.');grid;
legend('解的变化','种群均值的变化')
xlabel('迭代次数');ylabel('第一目标函数值');
figure(2);clf;
plot(trace2(:,1));hold on;
plot(trace2(:,2),'-.');
plot(trace2(:,1),'.');
plot(trace2(:,2),'.');grid;
legend('解的变化','种群均值的变化')
xlabel('迭代次数');ylabel('第二目标函数值');
figure(3);clf;
plot(trace3(:,1));hold on;
plot(trace3(:,2),'-.');
plot(trace3(:,1),'.');
plot(trace3(:,2),'.');grid;
legend('解的变化','种群均值的变化')
xlabel('迭代次数');ylabel('目标函数值之和');
figure(4);clf;plot(f1(v));hold on;
plot(f2(v),'r-.');grid;