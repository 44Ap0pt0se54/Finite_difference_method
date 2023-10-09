close all
clear all

% Erreur en temps

ea = zeros(1,1800);

phi = ones(1800,2);
y = zeros(1800,1);

mi = [200:1:2000];

for M = 200 : 2000

    N=floor(sqrt(M/2)); % condition de convergence

    Uexp = explicit(M,N);

    Ua = analytic(M,N);

    Ugap = abs(Uexp - Ua);

    ea(M-199) = mean(mean(Ugap));

    y(M-199) = log(mean(mean(Ugap)));
    phi(M-199,1) = log(M);


end

figure(1)

loglog(mi,ea);
xlabel("subdivision en temps (Nt)");
ylabel("erreur absolue moyenne en temps");
title("Erreur en temps : Schéma implicite");

teta = inv(phi'*phi)*phi'*y;
