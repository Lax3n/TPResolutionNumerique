clear


n = 100
h = 1/(n+1) 

u = 2*ones(n,1);
v = - ones(n-1,1)
A = diag(u)+diag(v,1)+diag(v,-1);
A = (1/h**2)*A 

b = rand(n,1)  

y=A\b


//////////////////////////
//
//
//  fonction preconditionnement SSOR: 
//  en entree: n\in N*, r\in R^n , A = D-E-F \in M_n, w \in ]0,2[ 
//
//  On suppose Aii non nuls pour tous i=1,...,n
//
//  en sortie p tel que (D/w - E)D^-1(D/w - F)p = r 
//
function p=SSOR(n,A,w,r)
 
    p= zeros(n,1)
    z = zeros(n,1)    
    
    
    
    // Descente (D/w-E)z = r
 
    z(1) = w*r(1)/A(1,1)
    for i=2:n
        z(i)= w*(r(i)-A(i,1:i-1)*z(1:i-1))/A(i,i);

    end
     
    
    // Remontee (D/w-F)p = Dz
    
    p(n) =w*z(n) 
    for i=n-1:-1:1
        p(i) = w*(A(i,i)*z(i)-A(i,i+1:n)*p(i+1:n))/A(i,i)
    end
    
endfunction


// Richardson a pas variable preconditionné 

w= 1.5
eps=10.0E-10

x = zeros(n,1);
r = b-A*x;
nr0 = norm(r);
nr = nr0;
k = 0;


while (nr/nr0>eps)&(k<10000) then 
        
    k = k + 1;    
    res(k) = log(nr/nr0)/log(10);
    p=SSOR(n,A,w,r)
    q=A*p
    alpha=(p'*r)/(q'*p)
    x=x+alpha*p
    r=r-alpha*q
    nr=norm(r);
  
end

plot(res,'-xr')