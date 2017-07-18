function h = restrictHFieldToOmega(f)
f(f<1) = 1;
f(f>max) = max;
