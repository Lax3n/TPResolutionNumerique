clear;

function [x, nit, Fx] = newton(x0, F, JF, eps, nitmax)
// initialization
x = x0;
r = F(x);
nr = norm(r);
nr0= nr;
nit=0;

while( nr/nr0 > eps & nit < nitmax )

x = x -JF(x)\r;
r = F(x);
nr = norm(r);
nit = nit+1;

end
Fx = r;

endfunction