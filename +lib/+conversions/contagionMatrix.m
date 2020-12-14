function C = contagionMatrix(source)
%CONTAGIONMATRIX Gives the contagious matrix from the contact matrix
%                as defined in the source file.

M = lib.conversions.contactMatrix(source, [
    0  9
    10 19
    20 29
    30 39
    40 49
    50 59
    60 69
    70 inf
]);

C = M * diag(lib.loaders.static_susceptibility);

end

