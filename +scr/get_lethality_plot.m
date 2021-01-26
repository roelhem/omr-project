Lethality = [0.00022274,0.00022274,0.00022274,0.00022274,0.00022274,0.00221398,0.01351613,0.134212];
LethalityPercentages = Lethality * 100;

GroupCategory = P.InitialState.GroupCategory;

figure

b = bar(GroupCategory, LethalityPercentages);
b.FaceColor = 'flat';
b.CData = jet(length(Lethality));

title('Lethality per Age Group');
ylim([0 ceil(max(LethalityPercentages)) + 1]);
ytickformat('percentage');

m = length(Lethality);

for i = 1:m
    text(i, LethalityPercentages(i), ...
        strcat(num2str(LethalityPercentages(i), '%2.2f'), '%'), ...
        'vert', 'bottom', 'horiz', 'center' ...
    );
end