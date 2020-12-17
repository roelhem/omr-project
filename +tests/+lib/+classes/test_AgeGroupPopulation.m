%% Get data from cbs.
P = lib.classes.AgeGroupPopulation.fromCbs();
assert(P.size == 100);

%% Total when resize to one.
P = lib.classes.AgeGroupPopulation.fromCbs();
O = P.rescale([0 inf]);
V = ones(P.size, 1);
assert(O.weightedResize(V, P) == 100, strcat(...
    "Resized matrix is ",...
    num2str(O.weightedResize(V, P))...
));