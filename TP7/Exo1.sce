clear;

function   [x,nit,Fx,nrvec]=newton(x0,F,JF,eps,nitmax)
    //initialisation
    x=x0;
    r=F(x);
    nr=norm(r);
    nr0=nr;
    nrvec=[nr/nr0]
    nit=0;
    
    while (nr/nr0>eps & nit<nitmax)
        x=x-JF(x)\r;
        r= F(x);
        nr=norm(r)
        nrvec=[nrvec nr/nr0];
        nit=nit+1
    end
    Fx=r;
endfunction



function y=F(x)
    N=length(x);
    y=zeros(N,1);
    //y(1)=(N+1)^2*sinh(x(1))+2*x(1)^2 -x(2)^2 -1;  
    
    //for i=2:N-1
    //    y(i)=(N+1)^2*sinh(x(i))+2*x(i)^2-x(i+1)^2-x(i-1)^2-i^2;
        
    //end
    //y(N)=(N+1)^2*sinh(x(N))+2*x(N)^2-x(N-1)^2-N^2
    
    y=(N+1)^2*sinh(x)+2*x.^2-[x(2:N);0].^2-[0;x(1:N-1)].^2-([1:N]').^2;
    
endfunction    



function y=JF(x)
    N=length(x);
    y=zeros(N,N);
    y=diag((N+1)^2*cosh(x)+4*x)-2*diag(x(1:N-1),-1)-2*diag(x(2:N),+1);
    
endfunction




N = 5 ;
x0 = zeros(N,1) ;
eps = 1e-8 ;
nitmax = 50 ;
[x,nit,Fx,nrvec] = newton( x0, F, JF, eps, nitmax )
printf( 'La methode a converge en %d iterations, ||Fx|| = %e \n', nit, norm(Fx) )

scf(0) ; 
clf ;
plot( [0:nit], nrvec, 'r-+' )
xa = gca() ;
xa. log_flags = 'nln' ;
disp(Fx)
disp(x)