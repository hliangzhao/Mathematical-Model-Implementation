%% 接种疫苗函数，这是和遗传算法唯一不同的函数，可以用它代替染色体的交叉操作。

%chromosomeGroup:染色体组
%bachterinChromosome:疫苗染色体，即最好的染色体。从这个染色体上取疫苗
%parameter:接种疫苗的参数，即用什么方法接种
%inoculateChromosome:接种疫苗后的染色体
function inoculateChromosome=immunity(chromosomeGroup,bacterinChromosome,parameter)
[chromosomeGroupSum,chromosomeLength]=size(chromosomeGroup);
[row,bacterinChromosomeLength]=size(bacterinChromosome);
%chromosomeGroupSum:染色体的条数；chromosomeLength：染色体的长度
switch parameter
    case 1%随机选择染色体进行接种
        for i=1:chromosomeGroupSum
            %%%%%%%%%%%%从疫苗染色体上定位疫苗
            headDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上左边的点位
            if headDot==0%防止出现0点位
                headDot=1;
            end
            tailDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上右边的点位
            if tailDot==0%防止出现0点位
                tailDot=1;
            end
            if tailDot>headDot%防止右边的点位大于左边的点位
                dot=headDot;
                headDot=tailDot;
                tailDot=dot;
            end
            %%%%%%%%%%%%%接种
            randChromosomeSequence=round(rand(1)*chromosomeGroupSum);
            %随机产生1条染色体的序号，对这条染色体进行接种
            if randChromosomeSequence==0%防止产生0序号
                randChromosomeSequence=1;
            end
            inoculateChromosome(i,:)...%先把输入染色体传给输出
            =chromosomeGroup(randChromosomeSequence,:);
            %执行免疫，即从疫苗染色体上取出一段基因做疫苗，再注入到其它染色体中
            inoculateChromosome(i,headDot:tailDot)...
            =bacterinChromosome(1,headDot:tailDot);
        end
    case 2 %所有染色体挨个接种
        for i=1:chromosomeGroupSum
        %%%%%%%%%%%%从疫苗染色体上定位疫苗
            headDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上左边的点位
            if headDot==0%防止出现0点位
                headDot=1;
            end
            tailDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上右边的点位
            if tailDot==0%防止出现0点位
                tailDot=1;
            end
            if tailDot>headDot%防止右边的点位大于左边的点位
                dot=headDot;
                headDot=tailDot;
                tailDot=dot;
            end
            %%%%%%%%%%%%%接种
            inoculateChromosome(i,:)=chromosomeGroup(i,:);%先把输入染色体传给输出
            %执行免疫，即从疫苗染色体上取出一段基因做疫苗，再注入到其它染色体中
            inoculateChromosome(i,headDot:tailDot)...
            =bacterinChromosome(1,headDot:tailDot);
        end
    case 3 %接种位置是随机的
        for i=1:chromosomeGroupSum
        %%%%%%%%%%%%从疫苗染色体上定位疫苗
            headDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上左边的点位
            if headDot==0%防止出现0点位
                headDot=1;
            end
            tailDot=fix(rand(1)*bacterinChromosomeLength);
            %疫苗在染色体上右边的点位
            if tailDot==0%防止出现0点位
                tailDot=1;
            end
            if tailDot>headDot%防止右边的点位大于左边的点位
                dot=headDot;
                headDot=tailDot;
                tailDot=dot;
            end
            %%%%%%%%%%%%%在染色体上随机定位接种位置
            inoculateDot=fix(rand(1)*chromosomeLength);%随机选择染色体的接种点位
            if inoculateDot==0
                inoculateDot=1;
                inoculateChromosome(i,:)=chromosomeGroup(i,:);
                inoculateChromosome(i,inoculateDot:tailDot-headDot+1)...
                =bacterinChromosome(1,headDot:tailDot);
            elseif inoculateDot<=headDot
                inoculateChromosome(i,:)=chromosomeGroup(i,:);
                inoculateChromosome(i,inoculateDot:inoculateDot+tailDot-headDot)...
                =bacterinChromosome(1,headDot:tailDot);
            elseif (chromosomeLength-inoculateDot)>=(tailDot-headDot)
                inoculateChromosome(i,:)=chromosomeGroup(i,:);
                inoculateChromosome(i,inoculateDot:inoculateDot+tailDot-headDot)...
                =bacterinChromosome(1,headDot:tailDot);
            else
                inoculateChromosome(i,:)=chromosomeGroup(i,:);
                inoculateChromosome(i,headDot:tailDot)...
                =bacterinChromosome(1,headDot:tailDot);
            end
        end
end