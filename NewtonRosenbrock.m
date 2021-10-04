function [output] = NewtonRosenbrock(start, numiter)
temp = start;
xiterates = [temp(1)];
yiterates = [temp(2)];
for i = 1:numiter
    x = temp(1);
    y = temp(2);
    gradX = 400 * x^3 + 2*x - 2 - 400*x*y;
    gradY = -200*x^2 + 200*y;
    grad = [gradX gradY];
    hessian = [1200*x^2+2-400*y -400*x;-400*x 200];
    temp = temp - transpose(inv(hessian)*transpose(grad));
    xiterates = [xiterates;temp(1)];
    yiterates = [yiterates;temp(2)];
end
scatter(xiterates,yiterates, 'filled');
output = temp;


end