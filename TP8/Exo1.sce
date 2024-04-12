//
//  Modele stratigraphique de remplissage des bassins par les sediments 
//
//
//
//  calcul de la hauteur h(x,t) des sediments fonction de x et du temps 
//
//   


clear

//donnees 
Lx=2; // longueur du bassin en km 

//Discretisation en espace 
n=100;
dx=Lx/n;

//Discretisation en temps 
tf = 1.5; // Temps final 
ndt = 50
dt = tf/ndt; // Pas de temps 

// Newton 
eps = 1.0E-6; // critere d'arret du Newton 
kmax = 100; // nombre max d'iterations du Newton 


// coefficients de diffusion 
Km = 1;
Kc = 10;

// flux aux bords gauche et droit 
g0 = -20;
g1 = 0;

function s = f_hmer(t)
    s = 25 + 5*cos(12*t) ;
endfunction

function s = f_psi(u)
    if (u<0) then
        s = Kc*u;
    else 
        s = Km*u;
    end
endfunction

function s = f_psip(u)
    if (u<0) then
        s = Kc;
    else 
        s = Km;
    end
endfunction


function s = f_hinit(x)
    s = 25*exp(-8*x/Lx) + 10;
endfunction


function F = f(h,t)

for i=1:n
 b(i) = f_hmer(t)-h(i);   
 psi(i) = f_psi(b(i));
end
for i=2:n-1 
 F(i) = (2*psi(i)-psi(i+1)-psi(i-1))/dx**2  
end
F(1) = ( (psi(1)-psi(2))/dx - g0 )/dx
F(n) = ( (psi(n)-psi(n-1))/dx - g1 )/dx

endfunction


function A = Df(h,t)


A = zeros(n,n);
A = zeros(n,n);

for i=1:n
b(i) = f_hmer(t)-h(i);
psi(i) = f_psi(b(i));
dpsi(i) = - f_psip(b(i))
end
for i=2:n-1
A(i,i) = 2*dpsi(i)/dx**2
A(i,i+1) =- dpsi(i+1)/dx ** 2
A(i,i-1) = -dpsi(i-1)/dx ** 2
end
A(1,1) = dpsi(1)/dx ** 2
A(1,2) = -dpsi(2)/dx ** 2
A(n,n) = dpsi(n)/dx ** 2
A(n,n-1) = -dpsi(n-1)/dx ** 2

2*dpsi(i)/dx ** 2
       

endfunction





function [x,xn,yn,newt] = hnew(hold,t,dt,eps,kmax,newt)
///////////////////////////////////////////////////////////
//
//  Etant donnes hold, t, dt, 
//
// on cherche x tel que G(x) = 0 avec G(x) = x - hold - dt*f(x,t+dt)
//
//  en utilisant l'algorithme de Newton 
//
/////////////////////////////////////////////////////////    
    
// Algorithme de Newton pour la solution de G(x) = 0: 

x = hold; // initialisation avec x = hold 
r = x-hold-dt*f(x,t+dt);
nr = norm(r);
nr0 = nr;

k = 1;
xn(k) = newt;
yn(k) = 1;

while (nr/nr0>eps)&(k<kmax)&(nr>eps)
    k = k + 1;
    newt = newt + 1;
    xn(k) = newt;    

    A = eye(n,n) -dt*Df(x,t+dt);
    x =x -A\r;
    r = x-hold-dt*f(x,t+dt);
    nr = norm(r);
        

    
    yn(k) = nr/nr0;    
end

if (k==kmax) then
    'newton non converge',x,nr
    abort
end

    
endfunction



////////////////////////////////////////////////////////////////////////////////

for i=1:n
  X(i) = (i-1/2)*dx;
end  




//initialisation de h 
for i=1:n
 h(i)= f_hinit(X(i));    
end

hs(:,1) = h; // stockage de h init 

scf(2);
title('h')
plot (X,h,'r')  


//Boucle en temps  /////////////////////////////////////////////////////////////
t = 0;
newt = 1;
for nt=1:ndt


[h,xn,yn,newt] = hnew(h,t,dt,eps,kmax)


t = nt*dt 

 
hs(:,nt) = h; // on stocke la solution au temps n 
     
 //plot h 
 scf(2);
 title('h') 
 plot (X,h,'-r')  
 
 
//Newton convergence
scf(1);
epsm = 1E-14;
plot(xn,log10(yn+epsm),'r-*')
clear xn;
clear yn; 

end // fin boucle en temps /////////////////////////////////////////////////////




// on enleve les parties erodees 
for i=1:n
for m=2:ndt
  j = ndt-m+1;
  hs(i,j) = min(hs(i,j),hs(i,j+1));  
end
end

scf(4);
title('hs')
for nt=1:ndt
plot(X,hs(:,nt),'-b')
end




