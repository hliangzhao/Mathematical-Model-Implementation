clc, clear;
%% 图像分割问题
% 利用遗传算法进行图像分割的基本思想是：把图像中的像素按灰度值用域值M分成两类图像，
% 一类为目标图像 ，另一类为背景图像 。
% 图像 由灰度值在0～M之间的像素组成，图像 由灰度值在M+1～L-1（L为图像的灰度级数）之间的像素组成。
% 本例处理图像为256灰度级，将灰度分割域值编码为一个8位0、1二进制码串
load woman                      %调用MATLAB中Woman图像灰度值
figure(1);                      %画图
image(X);colormap(map);         
NIND=40;                        %个体数目(Number of individuals)
MAXGEN=50;                      %最大遗传代数(Maximum number of generations)
PRECI=8;                        %变量的二进制位数(Precision of variables)
GGAP=0.9;                       %代沟(Generation gap)
FieldD=[8;1;256;1;0;1;1];       %建立区域描述器(Build field descriptor)
Chrom=crtbp(NIND,PRECI);        %创建初始种群
gen=0;    
phen=bs2rv(Chrom,FieldD);       %初始种群十进制转换
ObjV=target(X,phen);            %计算种群适应度值
while gen<MAXGEN                %代沟(Generation gap)
    FitnV=ranking(-ObjV);       %分配适应度值(Assign fitness values)
    SelCh=select('sus',Chrom,FitnV,GGAP);     %选择
    SelCh=recombin('xovsp',SelCh,0.7);        %重组
    SelCh=mut(SelCh);                         %变异
    phenSel=bs2rv(SelCh,FieldD);              %子代十进制转换
    ObjVSel=target(X,phenSel);
    [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);  %重插入
    gen=gen+1;
end
[Y, I]=max(ObjV);
M=bs2rv(Chrom(I,:),FieldD);                   %估计域值
[m, n]=size(X);
for i=1:m
    for j=1:n
        if X(i,j)>M                           %灰度值大于域值时是白色
            X(i,j)=256;
        end
    end
end
figure(2)                                     %画出分割后目标图像
image(X);colormap(map);