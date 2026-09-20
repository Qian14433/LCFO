% DEMO: LCFO vs original CFO on a 30-D Sphere function (no CEC code needed)
%   For CEC2017 benchmarks, download the official competition code
%   (cec17_func.mexw64 + input_data) and call Get_CEC2017-style wrapper.
clear; clc;

dim = 30;
lb = -100; ub = 100;
MaxFEs = 30000;

fobj = @(x) sum(x.^2);          % Sphere, optimum 0 at origin

% --- our algorithm: LCFO (FE-driven, big-start pop0 = 100) ---
rng(1);
[best_pd, ~, curve_pd] = LCFO(100, MaxFEs, lb, ub, dim, fobj);

% --- original CFO (iteration-driven, T = (MaxFEs-pop)/(2*pop)) ---
pop = 30;
T = floor((MaxFEs - pop) / (2*pop));
rng(1);
[best_cfo, ~, curve_cfo] = CFO(pop, T, lb, ub, dim, fobj);

fprintf('Sphere D=30, budget %d FEs:\n', MaxFEs);
fprintf('  LCFO best error: %.4e\n', best_pd);
fprintf('  CFO       best error: %.4e\n', best_cfo);

figure;
semilogy(curve_pd, 'LineWidth', 1.5); hold on;
semilogy(linspace(1, numel(curve_pd), numel(curve_cfo)), curve_cfo, 'LineWidth', 1.5);
xlabel('generation'); ylabel('best objective (log)');
legend('LCFO', 'CFO (original)', 'Location', 'northeast');
title('Sphere function, D=30, MaxFEs=30000');
grid on;
