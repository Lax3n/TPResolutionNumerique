// Exercice 1
// 1
function s=F(y)
    s=y-exp(sin(y-%pi/4))
endfunction

// 2
y0=5
eps=10**-8
nitmax=50
k=0
yk=y0

while abs(F(yk))/abs(F(y0))>eps & k<nitmax
    yk=yk-((F(yk))**2)/(F(yk+F(yk))-F(yk))
    k=k+1
end

disp(yk)

// Exercice 2
n=10
x0=zeros(n,1)
r0=b-A.*x0
p0=r0
nr0=norm(r0)
k=0

xk=x0
rk=r0
pk=p0
nrk=nr0

while nrk/nr0>eps & k<kmax