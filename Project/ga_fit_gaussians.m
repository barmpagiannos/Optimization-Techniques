function [best, hist] = ga_fit_gaussians(U, y, M, bounds, ga, fit)
% GA to fit up to M Gaussians without using ga()

D = 5*M; % genes: [w, c1, c2, s1, s2] for each i

% --- init population
pop = init_population(ga.pop_size, M, bounds);

% --- evaluate
J = zeros(ga.pop_size,1);
for p = 1:ga.pop_size
    J(p) = fitness(pop(p), U, y, M, fit);
end

hist.best_f = zeros(ga.n_gen,1);

nElite = max(1, round(ga.elite_frac * ga.pop_size));

for g = 1:ga.n_gen
    % sort by fitness (minimization)
    [J, idx] = sort(J, 'ascend');
    pop = pop(idx);

    best = pop(1);
    hist.best_f(g) = J(1);

    % elitism
    newpop = pop(1:nElite);

    % generate rest
    while numel(newpop) < ga.pop_size
        % selection
        p1 = tournament_select(pop, J, ga.tourn_k);
        p2 = tournament_select(pop, J, ga.tourn_k);

        % crossover
        c1 = p1; c2 = p2;
        if rand < ga.p_cx
            [c1, c2] = blend_crossover(p1, p2, bounds);
        end

        % mutation
        if rand < ga.p_mut
            c1 = mutate_individual(c1, bounds, ga.mut_sigma);
        end
        if rand < ga.p_mut
            c2 = mutate_individual(c2, bounds, ga.mut_sigma);
        end

        newpop(end+1) = c1; %#ok<AGROW>
        if numel(newpop) < ga.pop_size
            newpop(end+1) = c2; %#ok<AGROW>
        end
    end

    pop = newpop(:);

    % evaluate new population
    for p = 1:ga.pop_size
        J(p) = fitness(pop(p), U, y, M, fit);
    end
end

% final best
[J, idx] = sort(J, 'ascend');
pop = pop(idx);
best = pop(1);

end

%% ---------- helpers ----------
function pop = init_population(N, M, bounds)
pop(N,1) = make_empty_individual(M);

for i = 1:N
    ind = make_empty_individual(M);
    ind.w  = bounds.w(1)  + (bounds.w(2) - bounds.w(1))*rand(M,1);
    ind.c1 = bounds.u1(1) + (bounds.u1(2) - bounds.u1(1))*rand(M,1);
    ind.c2 = bounds.u2(1) + (bounds.u2(2) - bounds.u2(1))*rand(M,1);
    ind.s1 = bounds.s1(1) + (bounds.s1(2) - bounds.s1(1))*rand(M,1);
    ind.s2 = bounds.s2(1) + (bounds.s2(2) - bounds.s2(1))*rand(M,1);

    pop(i) = ind;
end
end

function ind = make_empty_individual(M)
ind.w  = zeros(M,1);
ind.c1 = zeros(M,1);
ind.c2 = zeros(M,1);
ind.s1 = ones(M,1);
ind.s2 = ones(M,1);
end

function J = fitness(ind, U, y, M, fit)
yhat = model_predict(ind, U, M);
mse = mean((yhat - y).^2);

active = sum(abs(ind.w) > fit.tau_w);
J = mse + fit.lambda * active;
end

function parent = tournament_select(pop, J, k)
n = numel(pop);
idx = randi(n, k, 1);
[~, bestLocal] = min(J(idx));
parent = pop(idx(bestLocal));
end

function [c1, c2] = blend_crossover(p1, p2, bounds)
% BLX-alpha style crossover on each gene vector
alpha = 0.35;

[c1, c2] = deal(p1);

fields = {'w','c1','c2','s1','s2'};
for f = 1:numel(fields)
    a = p1.(fields{f});
    b = p2.(fields{f});
    lo = min(a,b); hi = max(a,b);
    range = hi - lo;

    childA = lo - alpha*range + (1+2*alpha)*range.*rand(size(a));
    childB = lo - alpha*range + (1+2*alpha)*range.*rand(size(a));

    c1.(fields{f}) = childA;
    c2.(fields{f}) = childB;
end

% clamp to bounds
c1 = clamp_individual(c1, bounds);
c2 = clamp_individual(c2, bounds);
end

function ind = mutate_individual(ind, bounds, mut_sigma)
% Gaussian perturbation (relative to variable ranges)
% Also occasionally "kill" a gaussian by shrinking weight.
M = numel(ind.w);

% mutate each gaussian with some probability
p_gene = 0.25;

for i = 1:M
    if rand < p_gene
        ind.w(i)  = ind.w(i)  + mut_sigma * (bounds.w(2)-bounds.w(1))  * randn;
        ind.c1(i) = ind.c1(i) + mut_sigma * (bounds.u1(2)-bounds.u1(1))* randn;
        ind.c2(i) = ind.c2(i) + mut_sigma * (bounds.u2(2)-bounds.u2(1))* randn;
        ind.s1(i) = ind.s1(i) + mut_sigma * (bounds.s1(2)-bounds.s1(1))* randn;
        ind.s2(i) = ind.s2(i) + mut_sigma * (bounds.s2(2)-bounds.s2(1))* randn;

        % occasional sparsity push: shrink weight
        if rand < 0.10
            ind.w(i) = 0.2 * ind.w(i);
        end
    end
end

% clamp
ind = clamp_individual(ind, bounds);
end

function ind = clamp_individual(ind, bounds)
ind.w  = min(max(ind.w,  bounds.w(1)),  bounds.w(2));
ind.c1 = min(max(ind.c1, bounds.u1(1)), bounds.u1(2));
ind.c2 = min(max(ind.c2, bounds.u2(1)), bounds.u2(2));
ind.s1 = min(max(ind.s1, bounds.s1(1)), bounds.s1(2));
ind.s2 = min(max(ind.s2, bounds.s2(1)), bounds.s2(2));
end
