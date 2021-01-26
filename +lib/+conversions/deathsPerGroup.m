function [R,T] = deathsPerGroup(varargin)
%CASESPERGROUP Gets the amount of corona cases per group in a certain
%period.

global rivm_cases;
global cbs_AgeGroupPopulation;

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

% Only show the deaths
T = T(T.Deceased == 'Yes', :);

% Get the groups.
S = groupsummary(T, 'Agegroup');

% Distribute the deaths below 50 over the other age groups.
belowFifty = S{S.Agegroup == '<50', 'GroupCount'};
S(S.Agegroup == '<50', :) = [];

targetScale = cbs_AgeGroupPopulation.rescale([
     0  9
    10 19
    20 29
    30 39
    40 49
]);

R(6:10,:) = S;
R{1:5,1} = targetScale.GroupCategory;
R{1:5,2} = targetScale.Population / sum(targetScale.Population) * belowFifty;


end