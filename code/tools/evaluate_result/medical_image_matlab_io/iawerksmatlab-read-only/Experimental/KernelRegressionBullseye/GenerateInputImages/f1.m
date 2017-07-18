function y = f1(x)
c1 = 0.5;
c2 = 0.8;
s  = -0.1;
y1 = 1.0;
y2 = 0.1;

A = [c1^3 c1^2 c1 1; c2^3 c2^2 c2 1; 3*c1^2 2*c1 1 0; 3*c2^2 2*c2 1 0];
b = [y1+s*c1 y2+s*c2 s s]';
p = A\b;

y = zeros(size(x));

y(x<=c1)=y1+s*x(x<=c1);

y(x>c1)  = p(1)*x(x>c1).^3 + p(2)*x(x>c1).^2 + p(3)*x(x>c1) + p(4); 

y(x>=c2) = y2+s*x(x>=c2);
