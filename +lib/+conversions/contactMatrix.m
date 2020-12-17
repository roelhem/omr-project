function M = contactMatrix(filename, ageGroups)
%CONTACTMATRIX Returns a contact matrix rescaled to the provided groups.

global cbs_AgeGroupPopulation;

% Get the contact matrix from the data file and it's scale.
D = lib.loaders.file_contactMatrix(filename);
m = height(D);

% Get the scales of both the data and the wanted matrix.
DScale = cbs_AgeGroupPopulation.rescale(D.Group.Boundaries);
MScale = cbs_AgeGroupPopulation.rescale(ageGroups);

% Calculate the contact matrix on the wanted scale.
A = table2array(D(1:m,1:m));
M = MScale.avgResize(MScale.weightedResize(A, DScale)', DScale);

end

