clear;


function s=f(x)
    s(1,1)=-5.*x(1)+2.*sin(x(1))+2.*cos(x(2))
    s(2,1)=-5.*x(2)+2.*sin(x(2))+2.*cos(x(1))
endfunction

disp(f([1,2]))


function s=Df(x)
    s(1,1)=-5+2.*cos(x(1))
    s(1,2)=-2.*sin(x(2))
    s(2,1)=-2.*sin(x(1))
    s(2,2)=-5+2.*cos(x(2))
endfunction

function s=g(x)
    s(1,1)=2.*sin(x(1))+2.*cos(x(2))
    s(2,1)=2.*sin(x(2))+2.*cos(x(1))
    s=1/5.*s
endfunction

epsilon=10**-15
// x1=zeros(2,1)
x1=[1000000;-10000]
r1=g(x1)-x1
nr1=norm(r1)
k=1
vect1=[]

function s=algo1()
    nrk=nr1
    xk=x1
    rk=r1
    while (nrk/nr1)>epsilon
        xk=g(xk)
        rk=g(xk)-xk
        nrk=norm(rk)
        k=k+1
        vect1=[vect1,log(nrk)]
    end
    disp(k)
    
    plot(1:length(vect1),vect1,'-r')
    s=xk

endfunction

disp(Df([1,2]))



r1=f(x1)
nr1=norm(r1)
k=1
vect2=[]

function s=algo2()
    nrk=nr1
    xk=x1
    rk=r1
    while nrk/nr1>epsilon
        xk=xk-Df(xk)\rk
        rk=f(xk)
        nrk=norm(rk)
        k=k+1
        vect2=[vect2,log(nrk)]
    end
    disp(k)
    
    plot(1:length(vect2),vect2,'-g')
    s=xk
endfunction

disp(algo1())
disp(algo2())