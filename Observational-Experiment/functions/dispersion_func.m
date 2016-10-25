function num=dispersion_func(k,sig,h)
g=9.81; %gravity
num=abs(sig^2-g*k*tanh(k*h));