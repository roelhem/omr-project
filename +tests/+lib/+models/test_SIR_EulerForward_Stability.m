%function [S,I,R] = SIR_EulerForward(Szero,Izero,Rzero,beta,alpha,dt,n)

%% Test model with time-step size 1 and compare to time-step size 2 (one group)
[S1,I1,R1] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 1, 100);
[S2,I2,R2] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 2, 50);

Sdiff = S1(1:2:end) - S2;
Idiff = I1(1:2:end) - I2;
Rdiff = R1(1:2:end)- R2;

Sreldiff = Sdiff ./ S2;
Ireldiff = Idiff ./ I2;
Rreldiff = Rdiff ./ R2;

plot(Sdiff)
hold on;
plot(Idiff)
plot(Rdiff)
hold off;
%% Again same thing but now time-step 0.5 vs. 1 (one group)
[S1,I1,R1] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 0.5, 100);
[S2,I2,R2] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 1, 50);

Sdiff = S1(1:2:end) - S2;
Idiff = I1(1:2:end) - I2;
Rdiff = R1(1:2:end)- R2;

Sreldiff = Sdiff ./ S2;
Ireldiff = Idiff ./ I2;
Rreldiff = Rdiff ./ R2;

plot(Sdiff)
hold on;
plot(Idiff)
plot(Rdiff)
hold off;
%% Again same thing but now time-step 0.1 vs. 0.22 (one group)
[S1,I1,R1] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 0.1, 100);
[S2,I2,R2] = lib.models.SIR_EulerForward(100, 1, 0, 1, 0.05, 0.2, 50);

Sdiff = S1(1:2:end) - S2;
Idiff = I1(1:2:end) - I2;
Rdiff = R1(1:2:end)- R2;

Sreldiff = Sdiff ./ S2;
Ireldiff = Idiff ./ I2;
Rreldiff = Rdiff ./ R2;

plot(Sdiff)
hold on;
plot(Idiff)
plot(Rdiff)
hold off;