%% Test model without transmissions (one group)
[S,I,R] = lib.models.SIR_EulerForward(10, 20, 30, 0, 0, 1, 100);
assert(isequal(S, ones(1,100)*10), 'S is not 10 everywhere.');
assert(isequal(I, ones(1,100)*20), 'I is not 20 everywhere.');
assert(isequal(R, ones(1,100)*30), 'R is not 30 everywhere.');

%% Test model without recovery (one group)
[S,I,R] = lib.models.SIR_EulerForward(10, 20, 0, 1, 0, 1, 100);
assert(not(isequal(S, ones(1,100)*10)), 'S is constant.');
assert(not(isequal(I, ones(1,100)*20)), 'I is constant.');
assert(isequal(R, zeros(1,100)), 'R is not zero.');

%% Test model without infected (one group)
[S,I,R] = lib.models.SIR_EulerForward(10, 0, 30, 1, 1, 1, 100);
assert(isequal(S, ones(1,100)*10), 'S is not 10 everywhere.');
assert(isequal(I, zeros(1,100)), 'I is constant.');
assert(isequal(R, ones(1,100)*30), 'R is not 30 everywhere.');

%% Independent groups give the same result when computed seperately.
[S,I,R] = lib.models.SIR_EulerForward(...
    [10,20],...     Initial values for S
    [1,5],...       Initial values for I
    [8,12],...      Initial values for R
    [1,0;0,0.8],... beta
    [0.4,1],...     alpha
    1,100 ...       dt and n
);
[S1,I1,R1] = lib.models.SIR_EulerForward(10,1,8, 1  ,0.4,1,100);
[S2,I2,R2] = lib.models.SIR_EulerForward(20,5,12,0.8,1  ,1,100);
assert(isequal(S(1,:), S1), 'S(1,:) differs from S1.');
assert(isequal(S(2,:), S2), 'S(2,:) differs from S2.');
assert(isequal(I(1,:), I1), 'I(1,:) differs from I1.');
assert(isequal(I(2,:), I2), 'I(2,:) differs from I2.');
assert(isequal(R(1,:), R1), 'R(1,:) differs from R1.');
assert(isequal(R(2,:), R2), 'R(2,:) differs from R2.');