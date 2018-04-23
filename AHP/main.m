%% 层次分析法
clc, clear

fid = fopen('data.txt', 'r');
n1 = input('请输入准则层指标个数:');
n2 = input('请输入方案层指标个数:');
%n1 = 6;  %准则层指标个数
%n2 = 3;  %方案层指标个数
a = [];
% 读取准则层判断矩阵
for i = 1:n1
    tmp = str2num(fgetl(fid));
    a = [a;tmp];
end
% 读取方案层判断矩阵：判断矩阵存入bi中
for i = 1:n1
    str1 = char(['b', int2str(i), '=[];']);
    str2 = char(['b', int2str(i), '=[b', int2str(i), ';tmp];']);
    eval(str1);
    for j = 1:n2
        tmp = str2num(fgetl(fid));
        eval(str2);
    end
end

% 求最大特征值及对应的归一化特征向量
[max(1),wA] = ahp(a);
for i = 1:n1
    str = char(['[max(', int2str(i+1), '),wb', int2str(i), ']=ahp(', 'b', int2str(i), ')']);
    eval(str);
end

%计算平均一致性指标
RIT=CalculationRI();

[RIA,CIA]=sglsortexamine(max(1),a,RIT);
for i = 1:n1
    str = char(['[RIb', int2str(i), ',CIb', int2str(i), ...
    ']=sglsortexamine(max(', int2str(i+1), '),b', int2str(i), ',RIT);']);
    eval(str);
end
%{
[RIb1,CIb1]=sglsortexamine(max(2),b1,RIT);
[RIb2,CIb2]=sglsortexamine(max(3),b2,RIT);
[RIb3,CIb3]=sglsortexamine(max(4),b3,RIT);
[RIb4,CIb4]=sglsortexamine(max(5),b4,RIT);
[RIb5,CIb5]=sglsortexamine(max(6),b5,RIT);
[RIb6,CIb6]=sglsortexamine(max(7),b6,RIT);
%}

dw=zeros(n2,n1);
% 写入循环
for i = 1:n1
    str = char(['dw(1:', int2str(n2), ',', int2str(i), ')=wb', int2str(i), ';']);
    eval(str);
end
%{
dw(1:3,1)=wb1;
dw(1:3,2)=wb2;
dw(1:3,3)=wb3;
dw(1:3,4)=wb4;
dw(1:3,5)=wb5;
dw(1:3,6)=wb6;
%}

CIC = []; RIC = [];
for i = 1:n1
    str1 = char(['CIC=[CIC, CIb', int2str(i), '];']);
    eval(str1);
    str2 = char(['RIC=[RIC, RIb', int2str(i), '];']);
    eval(str2);
end
%{
CIC=[CIb1;CIb2;CIb3;CIb4;CIb5;CIb6];
RIC=[RIb1;RIb2;RIb3;RIb4;RIb5;RIb6];
%}
tw=tolsortvec(wA,dw,CIC,RIC)';
disp('输出准则层对目标层的权重:');
disp(wA);
disp('输出准则层对方案层的权重:');
disp(dw);
disp('输出总排序权值:');
disp(tw);

res=tw';
[diffcult,num]=sort(res)
