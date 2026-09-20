function [Best_score, Best_pos, curve] = LCFO(pop, MaxFEs, lb, ub, dim, fobj)
% LCFO — Life-Cycle Caterpillar Fungus Optimizer (lifecycle CFO).
%   Models the post-infection stages of the O. sinensis lifecycle:
%     1) Spore differential dispersal  — phase-1 differential host attraction
%        (DE-inspired, honestly cited), partners from X and dormant archive;
%     2) Sclerotium dormancy archive   — replaced parents kept as dormant
%        "sclerotia", capacity fixed at pop0, DECOUPLED from population
%        reduction (survives the season, feeds difference vectors);
%     3) Colony colonization memory    — H=10 success-history slots each for
%        scale MF and crossover MCR (SHADE-style Lehmer-weighted updates);
%     4) Seasonal carrying capacity    — big-start linear population
%        reduction pop0 -> 4 over the budget "season";
%     5) Infection stage (CFO-native re-parasitism) — gated OFF: its measured
%        success collapses to ~0 after 10% of a limited budget (see paper
%        ablation: removing it is worth W20/T5/L4 at D100 vs the gated variant).
%   FE-driven. Call with pop = 100 (big start).

pop0 = pop;
X = initialization(pop, dim, ub, lb);
fitness = zeros(1, pop);
FEs = 0;
for i = 1:pop
    fitness(i) = fobj(X(i, :));
    FEs = FEs + 1;
end
[Best_score, best_idx] = min(fitness);
Best_pos = X(best_idx, :);

CR = 0.9;
Hmh = 10;  MFh = 0.5*ones(Hmh,1);  Mk = 1;   % P2: success-history slots
MCRh = 0.5*ones(Hmh,1);   % P3: CR memory slots
RA = 1;  RB = 1;            % EMA success counts (phase1, phase2)
A = zeros(0, dim);          % dormant-host archive
Amax = pop0;   % P1: fixed, never shrinks
maxT = ceil((MaxFEs - pop) / pop) + 2;
curve = zeros(1, maxT);

t = 0;
while FEs < MaxFEs
    t = t + 1;
    phi = FEs / MaxFEs;

    %% PR: big-start linear reduction
    Ntarget = max(4, round(4 + (pop0 - 4) * (1 - phi)));
    if size(X,1) > Ntarget
        [~, ord] = sort(fitness);
        keep = ord(1:Ntarget);
        X = X(keep,:);  fitness = fitness(keep);
    end
    pop = size(X,1);

    %% phase 1: differential host attraction with archive-fed difference
    succF = [];  succW = [];  succCR = [];
    [~, rankIdx] = sort(fitness);
    poolN = max(2, ceil(0.2 * pop));
    hostPool = rankIdx(1:poolN);
    for i = 1:pop
        Xh = X(hostPool(randi(poolN)), :);
        r1 = randi(pop);
        while r1 == i, r1 = randi(pop); end
        if ~isempty(A) && rand < size(A,1) / (pop + size(A,1))
            Xr2 = A(randi(size(A,1)), :);
        else
            r2 = randi(pop);
            while r2 == i || r2 == r1, r2 = randi(pop); end
            Xr2 = X(r2, :);
        end
        ri = randi(Hmh);                              % P2: draw slot
        CRi = min(max(MCRh(ri) + 0.1*randn, 0), 1);   % P3: CR from slot
        Fi = MFh(ri) + 0.1 * tan(pi * (rand - 0.5));  % Cauchy around slot
        while Fi <= 0
            Fi = MFh(ri) + 0.1 * tan(pi * (rand - 0.5));
        end
        Fi = min(Fi, 1);
        V = X(i,:) + Fi * (Xh - X(i,:)) + Fi * (X(r1,:) - Xr2);
        jrand = randi(dim);
        mask = rand(1, dim) < CRi;  mask(jrand) = true;
        X_new = X(i,:);  X_new(mask) = V(mask);
        X_new = min(max(X_new, lb), ub);
        fnew = fobj(X_new);
        FEs = FEs + 1;
        if fnew < fitness(i)
            A(end+1, :) = X(i,:);          % replaced parent becomes dormant
            if size(A,1) > Amax
                A(randi(size(A,1)), :) = [];
            end
            succF(end+1) = Fi;             %#ok<AGROW>
            succCR(end+1) = CRi;           %#ok<AGROW>
            succW(end+1) = fitness(i) - fnew;   %#ok<AGROW>
            X(i,:) = X_new;
            fitness(i) = fnew;
        end
    end

    %% AS: update scale memory from this generation's successes
    if ~isempty(succF)
        w = succW / (sum(succW) + eps);
        MFh(Mk) = sum(w .* succF.^2) / (sum(w .* succF) + eps);
        MFh(Mk) = min(max(MFh(Mk), 0.1), 1);
        if max(succCR) > 0
            MCRh(Mk) = sum(w .* succCR.^2) / (sum(w .* succCR) + eps);
        else
            MCRh(Mk) = 0;
        end
        Mk = mod(Mk, Hmh) + 1;
    end

    %% phase 2: CFO-native parasitism, gated by measured success rate
    pB = 0;   % LCFO: infection stage retired (see header note 5)
    succ2 = 0;
    [~, best_idx] = min(fitness);
    current_best_pos = X(best_idx, :);
    for i = 1:pop
        if rand >= pB
            continue;              % starved this generation
        end
        if rand < 0.5
            e = 3 * randn(1, dim);
            X_new = X(i, :) + e .* (rand * current_best_pos - rand * X(i, :));
        else
            f = rand * (1 - phi)^2;
            X_new = current_best_pos + f .* (rand * current_best_pos - rand * X(i, :));
        end
        X_new = min(max(X_new, lb), ub);
        fnew = fobj(X_new);
        FEs = FEs + 1;
        if fnew < fitness(i)
            fitness(i) = fnew;
            X(i, :) = X_new;
            succ2 = succ2 + 1;
        end
    end
    RA = 0.9 * RA + numel(succF);
    RB = 0.9 * RB + succ2;

    [current_best_fitness, current_best_idx] = min(fitness);
    if current_best_fitness < Best_score
        Best_score = current_best_fitness;
        Best_pos = X(current_best_idx, :);
    end
    curve(t) = Best_score;
end
curve = curve(1:t);
end
