% Plot result
dt = 1;

% Run the SIR_EulerForward model.
[S,I,R] = lib.models.SIR_EulerForward(...
    [2000 2000 2000],...                    Initial conditions for S.
    [   0    1    0],...                    Initial conditions for I.
    [   0    0    0],...                    Initial conditions for R.
    [
        0.2 0.4  0.1
        0.3 0.5  0.1
        0   0.01 0.1
    ],...                                 Transmision rates.
    [0.04 0.08 0.01],...                  Recovery rates.
    dt,...                                Timestep size.
    400 ...                               Sample-size
);

% Display the result.
fig = lib.plots.SIR_Lines(S,I,R,dt);
fig.Name = 'Simple SIR model.';