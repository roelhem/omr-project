%% Test vectorized run has same result. (without vaccination strategy, simple data)
InitialState = lib.classes.ModelState(['0-9';'10+'], [1000; 100]);
InitialState.Beta = [
    1 0.2
    0.1 2
];
InitialState.I = [0;1];

m = InitialState.m;
n = 100;
x = zeros(1,m*n);

vectResult = InitialState.vacc_vectorized_run(x);

MA = vectResult.getModelResult(1);
MB = InitialState.run(1,n);

assert(all(abs(MA.S - MB.S) < 0.001, 'all'));
assert(all(abs(MA.I - MB.I) < 0.001, 'all'));
assert(all(abs(MA.R - MB.R) < 0.001, 'all'));
assert(all(abs(MA.Nu - MB.Nu) < 0.001, 'all'));

assert(all(abs(MA.S - MB.S) < 0.001, 'all'));


%% Test vectorized run has same result. (without vaccination strategy)
InitialState = lib.classes.ModelState.standardData();

m = InitialState.m;
n = 100;
x = zeros(1,m*n);

vectResult = InitialState.vacc_vectorized_run(x);


MA = vectResult.getModelResult(1);
MB = InitialState.run(1,n);


assert(all(abs(MA.S - MB.S) < 0.001, 'all'));
assert(all(abs(MA.I - MB.I) < 0.001, 'all'));
assert(all(abs(MA.R - MB.R) < 0.001, 'all'));
assert(all(abs(MA.Nu - MB.Nu) < 0.001, 'all'));

assert(all(abs(MA.S - MB.S) < 0.001, 'all'));