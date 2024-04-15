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
yk=yO

while abs(F(yk))/abs(F(y0))>eps & k<nitmax
    yk=yk-((F(yk))**2)/(F(yk+F(yk))-F(yk))
    k=k+1
end