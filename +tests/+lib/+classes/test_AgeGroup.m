global cbs_populationTotal;

%% Test Resizing to one full group.
% Getting the total group.
ATotalGroup = lib.classes.AgeGroup("0+");
% Getting the group distribution from the cbs.
CBSGroups = lib.classes.AgeGroup(cbs_populationTotal.Ten.Group.Boundaries);
V = ones(CBSGroups.size, 1);
assert(ATotalGroup.resize(V, CBSGroups) == CBSGroups.size, "Resize doesn't give the the total amount of groups.");