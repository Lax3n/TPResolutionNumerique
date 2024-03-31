clear;
function s=f(x)
    s=exp(-cos(x))+x-1
endfunction

function s=fp(x)
    s=sin(x).*exp(-cos(x))+1
endfunction


x=-10:0.1:10
// scf(1)
// plot(x,f(x))
// scf(2)
// plot(x,fp(x),'-r')




function s=g(x)
    s=x-f(x)/fp(a)
endfunction

function s=gp(x)
    s=1-fp(x)/fp(a)
endfunction


// plot(x,g(x))


a=0.5692773 
epsilon=10**-14
x0=-3
x=x0
r0=g(x0)-x0
r=r0
nr0=norm(r0)
nrk=nr0
k=0
vect=[]


while (nrk/nr0)>epsilon & k<1000
    x=g(x)
    r=g(x)-x
    nrk=norm(r)
    k=k+1
    vect=[vect,log(nrk/nr0)]

end
disp(x)
disp(k)

plot(1:k,vect,'-r')
disp(gp(x))

// while (nrk/nr0)>epsilon & k<1000


