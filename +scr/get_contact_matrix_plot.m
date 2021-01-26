%% Get the original (source) contact matrix.
[~, original_C, original_groups] = lib.loaders.file_contactMatrix(CFile);
figure;
h = heatmap(original_C);
h.XLabel = "Questioned person";
h.YLabel = "Contacted person";
h.XDisplayLabels = original_groups;
h.YDisplayLabels = original_groups;
h.Title = "Original Contact matrix";

figure;
P.InitialState.plotHeatMap('Data', 'C');

figure;
P.InitialState.plotHeatMap('Data', 'Beta');