function [output] = camel(x,y)

% output = (4-(2.1*(x^2))+((x^4)/3))*(x^2) + (x*y) + (-4+4*(y^2))*(y^2);

% x1 = x;
% x2 = y;

% fact1 = sin(x1)*cos(x2);
% fact2 = exp(abs(1 - sqrt(x1^2+x2^2)/pi));

% output = -abs(fact1*fact2);

output = -0.5*sin(x-pi/2)*sin(y-pi/2)+0.03*((x^2 + y^2))^(0.5);

end