clear;


//--------------------------------------
mprintf( '\nEXO I\n' ) ;
//--------------------------------------
//EXO I


function s=F(y)
    s = y - exp(sin(y - %pi/4));
endfunction


y = 5 ;
eps = 1.e-8 ;
nitmax = 50 ;
r = F(y)
nr = abs(r)
nr1 = nr
k=0

while (nr/nr1>eps)&(k<nitmax);
  
    y = y - r^2/( F(y+r) - r )   
    r = F(y)
    nr = abs(r)
    k = k+1
    res1(k) = nr/nr1
end

mprintf( 'y = %f\n', y ) ;
mprintf( 'nit = %d\n', k ) ;
scf(1)
plot(log(res1)/log(10),'x-b')

//--------------------------------------




//--------------------------------------
mprintf( 'EXO II\n') ;
//--------------------------------------
//EXO II


n = 4;
A = [ 7, -2, -1, 0 ; -2, 8, -3, -1 ; -1, -3, 10, -2; 0,-1,-2,6 ] ;
b = [ -1 ; 3 ; 2; -4 ] ;

sol = A\b ; // solution de reference 

eps = 1.d-10
nitmax = 100 ;
x = zeros(n,1);
r = b-A*x;
p = r;
trr = r'*r
nr0 = norm(r);
nr = nr0;
k = 0;


while (nr/nr0>eps)&(k<nitmax) then 
        
  
    q = A*p;
    alpha = trr/(p'*q);
    x = x + alpha*p;
    r = r - alpha*q;
    nr = norm(r);
    trr1 = r'*r;
    beta1 = trr1/trr;
    trr = trr1;
    p = r + beta1*p;
    
    k = k + 1;    
    res2(k) = nr/nr0       
end



mprintf( 'x = [ %f, %f, %f, %f ]\n', x(1), x(2), x(3), x(4) ) ;
mprintf( 'nit = %d\n', k ) ;
scf(2)
plot(log(res2)/log(10),'x-b')

