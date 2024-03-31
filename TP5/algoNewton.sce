

clear;
function s=f(x)
    s=exp(-cos(x))+x-1
endfunction

function s=fp(x)
    s=sin(x).*exp(-cos(x))+1
endfunction


// a=1
epsilon=10**-14
x0=-3
x=x0
r0=f(x0)
r=r0
nr0=norm(r0)
nrk=nr0
k=0
vect=[]

while (nrk/nr0)>epsilon & k<1000
    x=x-r/fp(x)
    r=f(x)
    nrk=norm(r)
    k=k+1
    vect=[vect,log(nrk/nr0)]
end


disp(x)
disp(k)

plot(1:k,vect,'-g')