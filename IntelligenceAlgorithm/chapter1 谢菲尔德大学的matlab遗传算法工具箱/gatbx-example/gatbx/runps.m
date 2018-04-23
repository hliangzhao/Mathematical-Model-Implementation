function [x fval] =  runps(a,b,c,x0)
[x, fval] = patternsearch(@nestedfun,x0);
% Nested function that computes the objective function
    function y = nestedfun(x)
        y = (a - b*x(1)^2 + x(1)^4/3)*x(1)^2 + x(1)*x(2) + ...
         (-c + c*x(2)^2)*x(2)^2;
    end
end