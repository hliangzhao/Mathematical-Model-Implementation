%% 免疫算法
%这个算法几乎与遗传算法一样，只是多用了一个免疫函数
%免疫算法是遗传算法的变体，它不用杂交，而是采用注入疫苗的方法。
%疫苗是优秀染色体中的一段基因，把疫苗接种到其它染色体中

%注意：标准遗传算法的一个重要概念是，染色体是可能解的2进制顺序号，由这个序号在可能解的集合(解空间)中找到可能解
%这是免疫算法的主程序，它需要调用的函数如下。
%接种疫苗函数：
%function inoculateChromosome=immunity(chromosomeGroup,bacterinChromosome,parameter)
%parameter:1,随机制取染色体接种。2，每个染色体都接种。3，每个染色体都接种，但接种的位置是随机的
%这个函数实现对染色体的疫苗接种

%由染色体(可能解的2进制)顺序号找到可能解：
%x=chromosome_x(fatherChromosomeGroup,oneDimensionSet,solutionSum);

%把解代入非线性方程组计算误差函数：functionError=nonLinearSumError1(x);

%判定程是否得解函数：[solution,isTrue]=isSolution(x,funtionError,solutionSumError);

%选择最优染色体函数：
%[bestChromosome,leastFunctionError]=best_worstChromosome(fatherChromosomeGroup,functionError);

%误差比较函数：从两个染色体中，选出误差较小的染色体
%[holdBestChromosome,holdLeastFunctionError]...
% =compareBestChromosome(holdBestChromosome,holdLeastFunctionError,...
% bestChromosome,leastFuntionError)
%为染色体定义概率函数，好的染色体概率高，坏染色体概率低
%p=chromosomeProbability(functionError);

%按概率选择染色体函数：
%slecteChromosomeGroup=selecteChromome(fatherChromosomeGroup,p);

%父代染色体杂交产生子代染色体函数
%sonChrmosomeGroup=crossChromosome(slecteChromosomeGroup,2);

%防止染色体超出解空间的函数
%chromosomeGroup=checkSequence(chromosomeGroup,solutionSum)

%变异函数
%fatherChromosomeGroup=varianceCh(sonChromosomeGroup,0.8,solutionN);

%通过实验有如下结果：
%1。染色体应当多一些
%2。通过概率选择染色体，在迭代早期会有效选出优秀的染色体，使解的误差迅速降低，
%但随着迭代的进行，概率选择也会导致某种染色体在基因池中迅速增加，使染色体趋同，
%这就减少了物种的多样性，反而难以逼近解
%3。不用概率选择，仅采用染色体杂交，采用保留优秀染色体，也可以得到解
%4。单纯免疫效果不好，杂交+免疫效果比较好

%%%%%%%%%%%%%%%%%%%%%%%%程序开始运行

clear,clc;%清理内存，清屏
circleN=200;%迭代次数
format long

%%%%%%%%%%%%%%%构造可能解的空间，确定染色体的个数、长度
solutionSum=4;leftBoundary=-10;rightBoundary=10;
distance=1;chromosomeSum=500;solutionSumError=0.1;
%solutionSum:非线性方程组的元数(待解变量的个数)；leftBoundary:可能解的左边界；
%rightBoundary:可能解的右边界；distance:可能解的间隔，也是解的精度
%chromosomeSum:染色体的个数；solveSumError:解的误差
oneDimensionSet=leftBoundary:distance:rightBoundary;
%oneDimensionSet:可能解在一个数轴(维)上的集合
oneDimensionSetN=size(oneDimensionSet,2);%返回oneDimensionSet中的元素个数
solutionN=oneDimensionSetN^solutionSum;%解空间(解集合)中可能解的总数
binSolutionN=dec2bin(solutionN);%把可能解的总数转换成二进制数
chromosomeLength=size(binSolutionN,2);%由解空间中可能解的总数(二进制数)计算染色体的长度

%%%%%%%%%%%%%%%%程序初始化
%随机生成初始可能解的顺序号,+1是为了防止出现0顺序号
solutionSequence=fix(rand(chromosomeSum,1)*solutionN)+1;
for i=1:chromosomeSum%防止解的顺序号超出解的个数
if solutionSequence(i)>solutionN;
solutionSequence(i)=solutionN;
end
end
%染色体是解集合中的序号,它对应一个可能解
%把解的十进制序号转成二进制序号
fatherChromosomeGroup=dec2bin(solutionSequence,chromosomeLength);
holdLeastFunctionError=Inf;%可能解的最小误差的初值
holdBestChromosome=0;%对应最小误差的染色体的初值

%%%%%%%%%%%%%%%%%%开始计算
compute=1;
circle=0;
while compute%开始迭代求解
%%%%%%%%%%%%%1:由可能解的序号寻找解本身(关键步骤)
x=chromosome_x(fatherChromosomeGroup,oneDimensionSet,solutionSum);
%%%%%%%%%%%%%2：把解代入非线性方程计算误差
functionError=nonLinearSumError1(x);%把解代入方程计算误差
[solution,minError,isTrue]=isSolution(x,functionError,solutionSumError);
%isSolution函数根据误差functionError判定方程是否已经解开，isTrue=1,方程得解。solution是方程的解
    if isTrue==1
        '方程得解'
        solution
        minError
        return%结束程序
    end
    %%%%%%%%%%%%%3：选择最好解对应的最优染色体
    [bestChromosome,leastFunctionError]=best_worstChromosome(fatherChromosomeGroup,functionError);
    %%%%%%%%%%%%%4：保留每次迭代产生的最好的染色体
    %本次最好解与上次最好解进行比较，如果上次最好解优于本次最好解，保留上次最好解；
    %反之，保留本次最好解。保留的最好染色体放在holdBestChromosome中
    [holdBestChromosome,holdLeastFunctionError]...
    =compareBestChromosome(holdBestChromosome,holdLeastFunctionError,...
    bestChromosome,leastFunctionError);
    circle=circle+1
    %minError
    %solution
    holdLeastFunctionError
    if circle>circleN
        return
    end
    %%%%%%%%%%%%%%5:把保留的最好的染色体holdBestChromosome加入到染色体群中
    order=round(rand(1)*chromosomeSum);
    if order==0
        order=1;
    end
    fatherChromosomeGroup(order,:)=holdBestChromosome;
    functionError(order)=holdLeastFunctionError;
    
    %%%%%%%%%%%%%%%6:为每一条染色体(即可能解的序号)定义一个概率(关键步骤)
    %%%%%%%%%%%%%%%好的染色体概率高，坏的概率低。依据误差functionError计算概率
    [p,trueP]=chromosomeProbability(functionError);
    if trueP =='Fail'
        '可能解严重不适应方程，请重新开始'
        return%结束程序
    end
    %%%%%%%%%%%%%%%7：按照概率筛选染色体(关键步骤)
    %fa=bin2dec(fatherChromosomeGroup)%显示父染色体
    %从父染体中选择优秀染色体
    %selecteChromosomeGroup=selecteChromosome(fatherChromosomeGroup,p);
    %%%%%%%%%%%%%%%8：染色体杂交(关键步骤)
    %sle=bin2dec(selecteChromosomeGroup)%显示选择出来的解的序号(染色体)
    %用概率筛选出的染色体selecteChromosomeGroup进行杂交，产生子代染色体
    %sonChromosomeGroup=crossChromosome(selecteChromosomeGroup,2);
    %不用概率筛选出的染色体selecteChromosomeGroup进行杂交，而直接用上一代(父代)的
    sonChromosomeGroup=crossChromosome(fatherChromosomeGroup,2);
    %sonChromosomeGroup=immunity(fatherChromosomeGroup,holdBestChromosome,3);
    %把疫苗接种到其它染色体中
    sonChromosomeGroup=immunity(sonChromosomeGroup,holdBestChromosome,3);
    %cro=bin2dec(sonChromosomeGroup)%显示杂交后的子代染色体
    sonChromosomeGroup=checkSequence(sonChromosomeGroup,solutionN);%检查杂交后的染色体是否越界
    %%%%%%%%%%%%%%%9：变异
    %不杂交直接变异
    %fatherChromosomeGroup=varianceCh(fatherChromosomeGroup,0.1,solutionN);
    %杂交后变异
    fatherChromosomeGroup=varianceCh(sonChromosomeGroup,0.5,solutionN);
    fatherChromosomeGroup=checkSequence(fatherChromosomeGroup,solutionN);%检查变异后的染色体是否越界
end
