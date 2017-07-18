function v = regkernel(x,mu,s)

v = 1/(s*sqrt(2*pi)) * exp(-(x-mu).^2/(2*s^2));
