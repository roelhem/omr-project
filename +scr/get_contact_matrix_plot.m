cbs_AgeGroupPopulation;


D = datetime() - calmonths(1);
G = lib.classes.ModelState.STD_Boundaries;
CFile = 'Contact_matrix.csv';
scale = cbs_AgeGroupPopulation.rescale(G);


MS = lib.classes.ModelState(scale.Boundaries, scale.Population);
MS.loadContactMatrix(CFile);

MS.ReprNum
MS.plotHeatMap('Data', 'C')