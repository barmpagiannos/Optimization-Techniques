clc; clear; close all;

%% Problem setup (from the statement)
u1_min = -1; u1_max = 2;
u2_min = -2; u2_max = 1;

% Ground-truth function ONLY for data generation & evaluation (as requested)
f_true = @(u1,u2) sin(u1 + u2) .* sin(u2.^2);

%% Data: train / test (different sets!)
N_train = 600;
N_test  = 600;

rng(1);
Utr = [u1_min + (u1_max-u1_min)*rand(N_train,1), ...
       u2_min + (u2_max-u2_min)*rand(N_train,1)];
ytr = f_true(Utr(:,1), Utr(:,2));

rng(2);
Ute = [u1_min + (u1_max-u1_min)*rand(N_test,1), ...
       u2_min + (u2_max-u2_min)*rand(N_test,1)];
yte = f_true(Ute(:,1), Ute(:,2));

%% Model parameters
M = 15; % max Gaussians

% Bounds for genes
bounds.u1 = [u1_min u1_max];
bounds.u2 = [u2_min u2_max];

% sigma bounds (choose reasonable positive range)
bounds.s1 = [0.05 1.5];
bounds.s2 = [0.05 1.5];

% weight bounds (can be wider if needed)
bounds.w  = [-3 3];

%% GA hyperparameters
ga.pop_size   = 80;
ga.n_gen      = 120;
ga.elite_frac = 0.08;      % keep best 8%
ga.tourn_k    = 3;         % tournament size
ga.p_cx       = 0.85;      % crossover probability
ga.p_mut      = 0.20;      % mutation probability per child
ga.mut_sigma  = 0.12;      % mutation strength (relative scale)
ga.seed       = 42;

% Complexity penalty
fit.lambda = 3e-3;
fit.tau_w  = 0.05;

%% Run GA
rng(ga.seed);

[best, hist] = ga_fit_gaussians(Utr, ytr, M, bounds, ga, fit);

%% Evaluate
yhat_tr = model_predict(best, Utr, M);
yhat_te = model_predict(best, Ute, M);

mse_tr = mean((yhat_tr - ytr).^2);
mse_te = mean((yhat_te - yte).^2);

active = sum(abs(best.w) > fit.tau_w);

fprintf('\n=== RESULTS ===\n');
fprintf('Active Gaussians (|w|>%.3f): %d / %d\n', fit.tau_w, active, M);
fprintf('Train MSE: %.6f\n', mse_tr);
fprintf('Test  MSE: %.6f\n', mse_te);

%% Plots: fitness, prediction quality, surfaces
figure; 
plot(hist.best_f, 'LineWidth', 1.5); grid on;
xlabel('Generation'); ylabel('Best fitness J');
title('GA Convergence (best fitness)');

figure;
scatter(yte, yhat_te, 18, 'filled'); grid on;
xlabel('y true (test)'); ylabel('y hat (test)');
title(sprintf('Test scatter | MSE=%.5f | active=%d', mse_te, active));
refline(1,0);

% Surface comparison on grid (optional)
nGrid = 60;
[u1g,u2g] = meshgrid(linspace(u1_min,u1_max,nGrid), linspace(u2_min,u2_max,nGrid));
Ytrue = f_true(u1g,u2g);
Ug = [u1g(:) u2g(:)];
Yhat = reshape(model_predict(best, Ug, M), size(u1g));

figure;
surf(u1g,u2g,Ytrue); shading interp;
title('True f(u1,u2)'); xlabel('u1'); ylabel('u2'); zlabel('y');

figure;
surf(u1g,u2g,Yhat); shading interp;
title('Estimated \hat{f}(u1,u2) via GA'); xlabel('u1'); ylabel('u2'); zlabel('y');

%% Print final analytic expression (compact)
print_model_expression(best, M, fit.tau_w);
