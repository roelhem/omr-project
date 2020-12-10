function cells = ensureChar(cells, default)
%NULLTOVAL Summary of this function goes here
%   Detailed explanation goes here

notChars = cellfun(@(x)(not(ischar(x))), cells);
cells(notChars) = default;

end

