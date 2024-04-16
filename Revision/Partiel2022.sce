// Exercice 1
// 1.
u=1

for i=1:20
    u=1/(1+u**2)
end
disp(u)

// 2.
i=0
for j=1:30
    i=i+(j**2)
end
// 3.
x0=1
y0=3
z0=2
xk=x0
yk=y0
zk=z0

for g=1:20
    xk=3.*xk+8.*yk+4.*zk
    yk=5.*xk+16.*yk+8.*zk
    zk=6.*xk+4.*yk+7.*zk
end

// 4.

x0=1
y0=3
z0=2
xk=x0
yk=y0
zk=z0

for g=1:20
    xk2=3.*xk+8.*yk+4.*zk
    yk2=5.*xk+16.*yk+8.*zk
    zk2=6.*xk+4.*yk+7.*zk
    xk=xk2
    yk=yk2
    yk=zk2
end

// Exercice 2
clear;
// 1.
n=5
u=2.*eye(n,n)
v=-ones(n-1,1)
A=u+diag(v,1)+diag(v,-1)
disp(A)

// disp(A(3,2))

// 2.
function s=f(x)
    s=x.*(1-x)
endfunction

// i=1
// b=f(i/n+1)

// 3.

u=zeros(n,1)
N=10
eps=10**-8

for i=1:n
    b(i)=f(i/n+1)
end

function s=somme(i)
    for j=1:n
        s=s+A(i,j).*u(j)
    end
endfunction


k=0
while norm(b-A*u)<eps | k>=N
    for i=1:n
        u(i)=u(i)+(b(i)-somme(i))/A(i,i)
    k=k+1
    end
end

disp(u)

// 4.
// Il résoud des problème du style Ax=b dans Rn
// xk est donc le vecteur solution de cette equation


