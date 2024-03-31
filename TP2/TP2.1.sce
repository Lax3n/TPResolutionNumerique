// Initialisation des paramètres
b = 1;
a = 2;
alpha = 0.25;
N = 4;

// 1) Construction du tableau x avec les valeurs de la suite
x = zeros(N+1, 1); // Création d'un tableau de taille N+1 initialisé à zéro
for k = 0:N-1 do
    x(k+2) = x(k+1) + alpha * (b - a * x(k+1)); // Notez que Scilab commence l'indexation à 1
end

// 2) Affichage de la valeur de x(N+1) qui est x(N) dans notre boucle à cause de l'indexation à 1
disp(x(N+1));

// 3) Pour calculer x(N) sans stocker toutes les valeurs, on réutilise la variable x_scalar
x_scalar = 0;
for k = 0:N-1 do
    x_scalar = x_scalar + alpha * (b - a * x_scalar);
end
disp(x_scalar)

// 4) Calcul de l'erreur et représentation graphique
erreur = abs(x - b/a); // Erreur absolue de x(k) par rapport à b/a

// Affichage de l'erreur sur une échelle logarithmique
scf(); // Crée une nouvelle figure
semilogy(0:N, erreur);
xlabel('k');
ylabel('Erreur |x(k) - b/a|');
title('Erreur de la suite x(k) en échelle logarithmique');

L=1
n=4
dx=L/(n+1)

u= 2*ones(n,1)/(dx)
v=-ones(n-1,1)/(dx)
A= diag(u)+diag(v,1)+diag(v,-1)
disp(A)

// c=dx(1:1:n)
// disp(c)