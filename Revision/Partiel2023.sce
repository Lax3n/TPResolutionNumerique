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
n=4
eps=10**-10
kmax=100

A=[7, -2, -1, 0 ; -2 , 8, -3, -1; -1, -3, 10, -2; 0, -1, -2, 6]
b=[-1; 3; 2; -4]

x0=zeros(n,1)
r0=b-A*x0
p0=r0
nr0=norm(r0)
k=0

xk=x0
rk=r0
pk=p0
nrk=nr0

while nrk/nr0>eps & k<kmax 
    qk=A*pk
    ak=(rk'*rk)/(pk'*qk)
    xk=xk+ak.*pk
    rk2=rk
    rk=rk-ak.*qk
    bk=(rk'*rk)/(rk2'*rk2)
    pk=rk+bk.*pk
    nrk=norm(rk)
    k=k+1
end

disp(xk)
