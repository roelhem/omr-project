function fig = SIR_Lines(S,I,R,dt)
%SIR_LINES Plot the results of an SIR-system as lines.
%   This command will take as it's input an calculated SIR-system (from
%   for instance, the results of a model) and then returns a `figure` that
%   contains a 2d-plot of this system.

%% Determine the time axis.
n = width(S);
m = height(S);
T = dt*(0:n-1);

%% Setting up the figure.
fig = figure;
xlabel('time in days');
ylabel('amount of people');

colorsS = summer(m);
colorsI = copper(m);
colorsR = cool(m);

%% Plotting the lines in the figure.
hold on;
% Plotting the values
for i = 1:m
    % Creating the data-tips.
    dataTipRows = [
        dataTipTextRow(['Group ' num2str(i)], '')
        dataTipTextRow('S', S(i,:))
        dataTipTextRow('I', I(i,:))
        dataTipTextRow('R', R(i,:))
        dataTipTextRow('t', 'XData')
    ];
    
    % Plotting the S-line
    Sp = plot(T, S(i,:), '.-',...
        'Color', colorsS(i,:)...
    );
    Sp.DataTipTemplate.DataTipRows = dataTipRows;
    Sp.DataTipTemplate.DataTipRows(1).Label = ['S of "group ' num2str(i) '"'];
    
    % Plotting the I-line
    Ip = plot(T, I(i,:), '.-',...
        'Color', colorsI(i,:)...
    );
    Ip.DataTipTemplate.DataTipRows = dataTipRows;
    Ip.DataTipTemplate.DataTipRows(1).Label = ['I of "group ' num2str(i) '"'];
    
    % Plotting the R-line
    Rp = plot(T, R(i,:), '.-',...
        'Color', colorsR(i,:)...
    );
    Rp.DataTipTemplate.DataTipRows = dataTipRows;
    Rp.DataTipTemplate.DataTipRows(1).Label = ['R of "group ' num2str(i) '"'];
end
hold off;

end

