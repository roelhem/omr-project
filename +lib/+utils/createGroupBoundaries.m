function B = createGroupBoundaries(groupSize, cutoffAge)
%CREATEGROUPBOUNDARIES Creates group boundaries from a groupSize and a cutoffAge.
%   - Default groupSize = 10,
%   - Default cutoffAge = 70.

switch nargin
    case 2
        l = groupSize;
        c = cutoffAge;
    case 1
        l = groupSize;
        c = 70;
    otherwise
        l = 10;
        c = 70;
end

% Initialize the result
m = ceil(c / l) + 1;
B = zeros(m, 2);

% Fill B with the bounded groups.
for i = 1:(m-1)
    B(i,:) = [(i-1)* l, min(i*l, c) - 1];
end

% Fill B with the cutoff group.
B(m,:) = [c, inf];

end

