function [Uanaly] = analytic(M,N)

h = 1/(N-1);
k = 1/(M-1);

n = 2000;


xi=[0:h:1];

t=[0:k:1];

nk=[1:1:n];

w = nk*pi;

Kd = 2./(w.^3).*(((-1).^nk)-1);

D = diag(Kd);

Uanaly = sin(xi'*w)*D*exp(-(w.^2)'*t);
