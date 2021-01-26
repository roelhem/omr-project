function out = stratday(x, th_max)
%STRATDAY Summary of this function goes here
%   Detailed explanation goes here

if width(x) == 1
    x = x';
end
m = width(x);

vals = x / th_max;

groups_labels = repmat("Group ", m, 1);
for j = 1:m
    groups_labels(j) = strcat("[", num2str(j), "]: ", num2str(x(j)));
end

out = pie(vals, groups_labels);

end

