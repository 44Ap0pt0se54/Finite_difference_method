clear all
close all

%% Variables du problème

l_x = 1;
l_y = 1;

n_x = 100;
n_y = 100;

pas_x = l_x/(n_x-1); 
pas_y = l_y/(n_y-1);

a1 = 1/(pas_x*pas_x);
a2 = 1/(pas_y*pas_y);
a3 = -2*(a1 + a2);

xi=[0:pas_x:1];

%% Construction de la matrice A et F tel que AU = F

A = zeros(n_x*n_y);

F = zeros(n_x*n_y,1);

for i = 1: n_x

    for j = 1: n_y

        if i > 1 && i < n_x && j > 1 && j < n_y     % cas du noeud intérieur

            A(bijection(i,j,n_x), bijection(i-1,j,n_x)) = a1;
            A(bijection(i,j,n_x), bijection(i+1,j,n_x)) = a1;

            A(bijection(i,j,n_x), bijection(i,j,n_x)) = a3;
            
            A(bijection(i,j,n_x), bijection(i,j-1,n_x)) = a2;
            A(bijection(i,j,n_x), bijection(i,j+1,n_x)) = a2;

            F(bijection(i,j,n_x)) = f((i-1)*pas_x, (j-1)*pas_y, l_x, l_y);

        elseif i == 1 || i == n_x || j==n_y         % CI de type Dirichlet

            A(bijection(i,j,n_x), bijection(i,j,n_x)) = 1;

            F(bijection(i,j,n_x)) = 0;

        elseif j == 1 && i>1 && i<n_x               % CI de type Neumann

            F(bijection(i,j,n_x)) = -g((i-1)*pas_x,l_x,l_y);
            
            % Methode d'ordre 1
            
            % A(bijection(i,j,n_x),bijection(i,1,n_x)) = -1/pas_y;
            % A(bijection(i,j,n_x),bijection(i,2,n_x)) = 1/pas_y;

            % Methode d'ordre 2

            A(bijection(i,j,n_x),bijection(i,1,n_x)) = -3/(2*pas_y);
            A(bijection(i,j,n_x),bijection(i,2,n_x)) = 4/(2*pas_y);
            A(bijection(i,j,n_x),bijection(i,3,n_x)) = -1/(2*pas_y);


        end
    end
end

%% Construction de U

U = inv(A)*F;
U=reshape(U,n_x,n_y);

figure(11)
surf(xi,xi,U);
title("Solution de l'équation de Poisson dans le plan")
hold on

% U theorique

Uth = zeros(n_x,n_y);
for i = 1: n_x
    for j = 1: n_y

        Uth(i,j) = sin((pi*(i-1)*pas_x)/l_x)*sin((pi*(j-1)*pas_y)/l_y);

    end
end



% Erreur d'approximation Neumann ordre 1



eps = norm(Uth - U);





