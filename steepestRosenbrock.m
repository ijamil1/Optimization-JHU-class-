function [output] = steepestRosenbrock(start, numiter)
temp = start;
xiterates = [temp(1)];
yiterates = [temp(2)];
for i = 1:numiter
    x = temp(1);
    y = temp(2);
    gradX = 400*x^3 + 2*x - 2 - 400*x*y;
    gradY = -200*x^2 + 200*y;
    a = -1 * gradX;
    b = -1 * gradY;
    alpha = optimalAlphaSteepestDescentRosenbruckFunc(temp, [a b]);
    temp(1) = temp(1) + alpha*a;
    temp(2) = temp(2) + alpha*b;
    xiterates = [xiterates;temp(1)];
    yiterates = [yiterates;temp(2)];
end
scatter(xiterates,yiterates);
output = temp;


end