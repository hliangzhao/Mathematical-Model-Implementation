function state = gaplotchange(options, state, flag)
% GAPLOTCHANGE Plots the logarithmic change in the best score from the
% previous generation.绘图函数gaplotchange画出前一代个体最佳值变化的图形
%   
persistent last_best % Best score in the previous generation 前一代个体的最佳值

if(strcmp(flag,'init')) % Set up the plot 开始绘图
    set(gca,'xlim',[1,options.Generations],'Yscale','log');
    hold on;
    xlabel Generation
    title('Change in Best Fitness Value')
end

best = min(state.Score); % Best score in the current generation 当前代个体的最佳值
if state.Generation == 0 % Set last_best to best. 设置last_best为最佳
    last_best = best;
else
 change = last_best - best; % Change in best score 个体最佳值的改变
 last_best=best;
 plot(state.Generation, change, '.r');
 title(['Change in Best Fitness Value'])
end