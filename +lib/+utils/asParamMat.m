function out = asParamMat(x, m)
%ASPARAMMAT Converts a given input `x` to a square matrix of size m.
%   - If x is already a square matrix of size m, it will just return that
%   matrix. 
%   - If x is a vector of size m, it will convert it to a matrix that
%   has the values of this matrix as it's diagonal.
%   - If x is a scalar, it will return a diagonal matrix that has the
%   value x on it's diagonal.
%   - This method throws an error otherwise.

if width(x) == m && height(x) == m
    out = x;
elseif length(x) == m && min(size(x)) == 1
    out = diag(x);
elseif length(x) == 1
    out = eye(m)*x;
else
    error('x cannot be converted to a ParamMat of size m.');
end
end

