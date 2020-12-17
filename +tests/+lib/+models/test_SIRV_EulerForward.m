%function [S,I,R,V] = SIRV_EulerForward(Szero,Izero,Rzero,Vzero,beta,alpha,VaccinationFrac,dt,n)
%% Test model with 1% vaccination in each day
VaccinationFrac = ones(1,100) .* 0.01;
[S,I,R,V] = lib.models.SIRV_EulerForward(100, 1, 0, 0, 1, 1, VaccinationFrac, 1, 100);

hold on;
plot(S)
plot(I)
plot(R)
plot(V)
hold off;

%% Test model with 0.5% vaccination in each day
VaccinationFrac = ones(1,100) .* 0.005;
[S,I,R,V] = lib.models.SIRV_EulerForward(100, 1, 0, 0, 1, 1, VaccinationFrac, 1, 100);

hold on;
plot(S)
plot(I)
plot(R)
plot(V)
hold off;