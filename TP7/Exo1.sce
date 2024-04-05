clear;

function [x, nit, Fx,nrvec] = newton(x0, F, JF, eps, nitmax)
    // initialization
    x = x0;
    r = F(x);
    nr = norm(r);
    nr0= nr;
    nit=0;
    nrvec=[nr/nr0]
    while( nr/nr0 > eps & nit < nitmax )

        x = x -JF(x)\r;
        r = F(x);
        nr = norm(r);
        nrvec=[nrvec nr/nr0]
        nit = nit+1;

    end
    Fx = r;

endfunction


function y=F(x)
    n=length(x)
    y=zeros(n,1)
    // y(1)=(n+1)^2*sinh(x(1))+2*(x(1))^2-(x(2))^2-1
    // for i=2:n-1 
    //     y(i)=(n+1)^2*sinh(x(i))+2*(x(i))^2-(x(i+1))^2-(x(i-1))^2-i^2
    // end
    // y(n)=(n+1)^2*sinh(x(n))+2*(x(n))^2-(x(n-1))^2-n^2

    y=(n+1)^2*sinh(x) +2*x.^2-[x(2:n);0].^2 - [0;x(1:n-1)].^2 - ([1:n]').^2;

endfunction

// disp(F([1;2;3]))

function y=JF(x)
    N=length(x);
    y=zeros(N,N)
    y=diag((N+1)^2*cosh(x)+4*x) -2* diag(x(1:N-1),-1)-2* diag(2:N,+1)
endfunction


// disp(JF([1;2;3]))


N = 5 ;
x0 = zeros(N,1) ;
eps = 1e-8 ;
nitmax = 50 ;
[x,nit,Fx,nrvec] = newton( x0, F, JF, eps, nitmax )
printf( 'La methode a converge en %d iterations, ||Fx|| = %e \n', nit, norm(Fx) )
scf(0) ; clf ;
plot( [0:nit], nrvec, 'r-+' )
xa = gca() ;
xa. log_flags = 'nln' ;

disp(x)