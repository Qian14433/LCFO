function [Best_score, Best_pos, curve] = CFO(pop, Max_iter, lb, ub, dim, fobj)

X = initialization(pop, dim, ub, lb); 

fitness = zeros(1, pop);
for i = 1:pop
    fitness(i) = fobj(X(i, :));
end

[Best_score, first_best_idx] = min(fitness);
Best_pos = X(first_best_idx, :);

curve = zeros(1, Max_iter); 

for t = 1:Max_iter
    
    [~, best_idx] = min(fitness);
    current_best_pos = X(best_idx, :);
    [sorted_fitness, sort_indices] = sort(fitness);
    sorted_X = X(sort_indices, :);
    
    X_new_exploration = zeros(pop, dim); 
    if rand < 0.5
        r1 = rand;
        beta = 2 * cos(pi * r1) * abs((t / Max_iter)^(r1 * randi([1 2])));
        X_new_exploration(1, :) = current_best_pos - rand(1, dim) .* (current_best_pos - sorted_X(1, :)) + beta .* (current_best_pos - sorted_X(1, :));
    else
        r = rand(1, dim);
        alpha = 2.5 * r .* abs(cos(pi * r));
        X_new_exploration(1, :) = sorted_X(1, :) - rand(1, dim) .* (current_best_pos - sorted_X(1, :)) + alpha .* (current_best_pos - sorted_X(1, :));
    end
    for i = 2:pop
        if rand < 0.5
            r1 = rand;
            beta = 2 * cos(pi * r1) * abs((t / Max_iter)^(r1 * randi([1 2])));
            X_new_exploration(i, :) = current_best_pos - rand(1, dim) .* (sorted_X(i-1, :) - sorted_X(i, :)) + beta .* (current_best_pos - sorted_X(i, :));
        else
            r = rand(1, dim);
            alpha = 2.5 * r .* abs(cos(pi * r));
            X_new_exploration(i, :) = sorted_X(i, :) - rand(1, dim) .* (sorted_X(i-1, :) - sorted_X(i, :)) + alpha .* (current_best_pos - sorted_X(i, :));
        end
    end

    for i = 1:pop
        X_new_exploration(i, :) = min(max(X_new_exploration(i, :), lb), ub);
        new_fitness = fobj(X_new_exploration(i, :));
        
        if new_fitness < sorted_fitness(i)
            original_index = sort_indices(i);
            X(original_index, :) = X_new_exploration(i, :);
            fitness(original_index) = new_fitness;
        end
    end

    [~, best_idx] = min(fitness);
    current_best_pos = X(best_idx, :);
    
    for i = 1:pop
        if rand < 0.5
            e = 3 * randn(1, dim);
            X_new = X(i, :) + e .* (rand * current_best_pos - rand * X(i, :));
        else
            f = rand * ((1 / Max_iter^2) * t^2 - 2 / Max_iter * t + 1);
            X_new = current_best_pos + f .* (rand * current_best_pos - rand * X(i, :));
        end
        
        X_new = min(max(X_new, lb), ub);
        new_fitness = fobj(X_new);
        
        if new_fitness < fitness(i)
            fitness(i) = new_fitness;
            X(i, :) = X_new;
        end
    end
    
    [current_best_fitness, current_best_idx] = min(fitness);
    if current_best_fitness < Best_score
        Best_score = current_best_fitness;
        Best_pos = X(current_best_idx, :);
    end
    
    curve(t) = Best_score;
end
end

