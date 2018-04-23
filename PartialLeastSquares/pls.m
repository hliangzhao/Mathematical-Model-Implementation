%% 使用偏最小二乘回归分析数据
clc, clear

load data.txt     % 注意此处对数据的要求：前三列 为自变量，后三列为因变量
mu = mean(data);  % 均值
sig = std(data);  % 标准差

rr = corrcoef(data);       % 求解相关系数矩阵
std_data = zscore(data);   % 数据标准化
n = 3; m = 3;              % n为自变量的个数，m为因变量的个数
x0 = data(:,1:n);
y0 = data(:,n+1:end);
e0 = std_data(:,1:n);
f0 = std_data(:,n+1:end);
num = size(e0,1);         % 样本点的个数

chg = eye(n);              % w到w*变换矩阵的初始化
for i = 1:n                
    % 计算w，w*和t的分向量
    matrix = e0'*f0*f0'*e0;
    [vec,val] = eig(matrix);
    val = diag(val);
    [val,ind] = sort(val,'descend');
    w(:,i) = vec(:,ind(1));           % 提出最大特征值对应的特征向量
    w_star(:,i) = chg*w(:,i);       % 计算w*的取值
    t(:,i) = e0*w(:,i);             % 计算成分ti的得分
    alpha = e0'*t(:,i)/(t(:,i)'*t(:,i));
    chg = chg*(eye(n)-w(:,i)*alpha');
    e = e0-t(:,i)*alpha';
    e0 = e;

    % 计算ss(i)的值
    beta = [t(:,1:i),ones(num,1)]\f0;   % 回归方程的系数
    beta(end,:) = [];                       % 删除系数的常数项
    residual = f0-t(:,1:i)*beta;        % 求残差矩阵
    ss(i) = sum(sum(residual.^2));          % 求误差平方和

    % 计算press(i)的值
    for j = 1:num
        t1 = t(:,1:i); f1 = f0;
        discard_t = t1(j,:); discard_f = f1(j,:);    % 保存舍去的第j个样本点
        t1(j,:) = []; f1(j,:) = [];                  % 删除第j个观测值
        beta1 = [t1, ones(num-1,1)]\f1;             % 求回归分析的系数
        beta1(end,:) = [];                            % 删除回归系数的常数项
        residual = discard_f-discard_t*beta1;        % 求残差向量
        press_i(j) = sum(residual.^2);
    end

    press(i) = sum(press_i);
    if i > 1
        Q_h2(i) = 1-press(i)/ss(i-1);
    else
        Q_h2(1) = 1;
    end

    if Q_h2(i) < 0.0975
        fprintf('主成分个数r=%d', i);
        r = i;
        break;
    end
end

beta_z = [t(:,1:r),ones(num,1)]\f0;              % 求y关于t的回归系数
beta_z(end,:) = [];                                  % 删除常数项
coeff = w_star(:,1:r)*beta_z;                      % 求y关于x的回归系数(是对标准化后的数据而言的)，每一列为一个回归方程

mu_x = mu(1:n);  mu_y = mu(n+1:end);
sig_x = sig(1:n);  sig_y = sig(n+1:end);

for i = 1:m
    ch0(i) = mu_y(i)-mu_x./sig_x*sig_y(i)*coeff(:,i); % 计算原始数据的回归方程系数的常数项
end

for i = 1:m
    coeff_origin(:,i) = coeff(:,i)./sig_x'*sig_y(i);  % 计算原始数据的回归方程系数，每一列为一个回归方程
end

sol = [ch0;coeff_origin];      % 回归方程系数

% 绘制预测图
ch0 = repmat(ch0, num, 1);
y_hat = ch0 + x0*coeff_origin;
y1_max = max(y_hat);
y2_max = max(y0);
y_max = max([y1_max;y2_max]);
residual = y_hat - y0;

subplot(2, 2, 1);
plot(0:y_max(1), 0:y_max(1), y_hat(:,1), y0(:,1), '*');

subplot(2, 2, 2);
plot(0:y_max(2), 0:y_max(2), y_hat(:,2), y0(:,2), 'O');

subplot(2, 2, 3);
plot(0:y_max(3), 0:y_max(3), y_hat(:,3), y0(:,3), 'H');