function [R,T] = casesPerGroup(varargin)
%CASESPERGROUP Gets the amount of corona cases per group in a certain
%period.

global rivm_cases;

% Filter the dataset.
T = rivm_cases;
if (nargin == 1)
    T = T(T.Date_statistics <= datetime(varargin{1}), :);
elseif (nargin >= 2)
    filter = and(...
        T.Date_statistics >= datetime(varargin{1}),...
        T.Date_statistics <= datetime(varargin{2})...
    );
    T = T(filter, :);
end

% Get the groups.
R = groupsummary(T, 'Agegroup');

% Remove the <50 and Unknown groups, as they are very small.
R = R(and(R.Agegroup ~= '<50', R.Agegroup ~= 'Unknown'), :);

end