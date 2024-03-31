n=10
b=ones(n,1)
a=2
alpha=0.25
L=1
dx=L/(n+1)
u= 2*ones(n,1)/(dx)
v=-ones(n-1,1)/(dx)
A= diag(u)+diag(v,1)+diag(v,-1)
X=zeros(n,1)
r=b-A*X
epsilon=1e-10
nr0=norm(r)
nr=nr0
k=0
vp=spec(A)
vpmin=min(vp)
vpmax=max(vp)
alphaopti=2/(vpmin+vpmax)

while((nr/nr0)>epsilon & k<10000) do
    k=k+1
    res(k)=log(nr/nr0)/log(10)
    p=A*r
    alphavariable=(r'*r)/(p'*r) //pas variable
    X=X+alphavariable*r
    r=r-alphavariable*p
    nr=norm(r)
    // res(k)=log(nr/nr0)/log(10)
    
end

// disp(X)
disp(k)
// disp(res)
plot(1:k,res,"-xr")