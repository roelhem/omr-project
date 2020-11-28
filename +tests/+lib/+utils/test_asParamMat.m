%% Test scalar input
y = lib.utils.asParamMat(2,3);
display(y);
assert(isequal(y, [2,0,0;0,2,0;0,0,2]), 'Unexpected result from scalar input')

%% Test vector input
y = lib.utils.asParamMat([1,2,3],3);
assert(isequal(y, [1,0,0;0,2,0;0,0,3]), 'Unexpected result from vector input')

%% Test matrix input
x = [1,2,3;4,5,6;7,8,9];
y = lib.utils.asParamMat(x,3);
assert(isequal(y, x), 'Unexpected result from matrix input')