%% Scalar linear growth: f(W) = W + 1
outputDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(outputDir));
addpath(fullfile(projectDir, 'dependencies'), '-begin');

model = struct('A',0, 'B',0, 'L',1, 'R',1, ...
               'P',0, 'tau',1, 'sigma',0);

cfg = struct('Ns',201, 'maxIter',2e5, 'maxWNorm',1e12, ...
             'relTol',1e-6, 'progressInterval',5, ...
             'odeSolver','ode15s', 'odeMaxTime',5000, ...
             'odeRelTol',1e-9, 'odeAbsTol',1e-11);

criterion = '||Wnew-Wold||_F <= relTol*||Wold||_F';
matlabVersion = version;
solverFile = fullfile('dependencies', 'riccati.m');

result = riccati(model, cfg);

save(fullfile(outputDir, 'results.mat'), ...
     'model', 'cfg', 'criterion', 'result', ...
     'matlabVersion', 'solverFile');

disp(result);
