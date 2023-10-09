clear all
close all

%% Variables du problème

mi = [10:50];

eps = zeros(1,40)

phi = ones(40,2);
y = zeros(40,1);

for Nxy = 10 : 50

    eps(Nxy-9) = norm(Uth(Nxy) - Uxy(Nxy));

    y(Nxy-9) = log(norm(Uth(Nxy) - Uxy(Nxy)));
    phi(Nxy-9,1) = log(Nxy);

end

figure(1)

loglog(mi,eps);
xlabel("subdivision x et y");
ylabel("erreur approximation");
title("Erreur d'approximation : Neumann ordre 2");

teta = inv(phi'*phi)*phi'*y;
    
