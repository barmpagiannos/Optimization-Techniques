
% Τεχνικές Βελτιστοποίησης
% 2η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc; clear; close all;

% Ορίζω τη συνάρτηση προς βελτιστοποίηση.
f = @(x) x(1)^3 * exp(-x(1)^2 - x(2)^4);

% Κλίση της συνάρτησης.
% x(1) = x, x(2) = y.
grad_f = @(x) [
    (3*x(1)^2 - 2*x(1)^4) * exp(-x(1)^2 - x(2)^4); % df/dx.
    -4 * x(1)^3 * x(2)^3 * exp(-x(1)^2 - x(2)^4) % df/dy.
];

% Τα αρχικά σημεία.
initial_points = [0, 0; -1, -1; 1, 1]';
step_strategies = {'constant', 'min', 'armijo'};
strategy_titles = {'Σταθερό Βήμα (0.1)', 'Ελαχιστοποίηση', 'Armijo Rule'};
epsilon = 1e-3; % ε = 0.001.
max_iter = 100; % μέγιστος αριθμός επαναλήψεων.

% Φτιάχνουμε το γράφημα των ισοϋψών καμπυλών.
[X, Y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
Z = (X.^3) .* exp(-X.^2 - Y.^4);

% Εκτέλεση της μεθόδου μέγιστης καθόδου για κάθε αρχικό σημείο.
for i = 1:3 % Για κάθε αρχικό σημείο (0,0), (-1,-1), (1,1).
    
    x0 = initial_points(:, i);
    point_str = sprintf('(%d, %d)', x0(1), x0(2));
    colors = ['r', 'g', 'b']; % Χρώματα.
    
    fprintf('\n--- Αρχικό Σημείο: %s ---\n', point_str);
    
    % Plot: πορεία σύγκλισης.
    figure('Name', ['Πορεία Σύγκλισης - Start: ' point_str], 'NumberTitle', 'off', 'Position', [100, 500, 1200, 350]);
    
    for j = 1:3 % Για κάθε διαφορετική επιλογή βήματος.
        step_type = step_strategies{j};
        subplot(1, 3, j);
        contourf(X, Y, Z, 20); hold on; 

        % Σχεδίαση αρχικού σημείου - ΑΣΤΕΡΙ.
        plot(x0(1), x0(2), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y'); 

        % Κλήση μεθόδου μέγιστης καθόδου.
        [x_hist, f_hist, k] = steepest_descent(f, grad_f, x0, step_type, epsilon, max_iter);
        plot(x_hist(1, :), x_hist(2, :), '.-', 'Color', colors(j), 'LineWidth', 1.5, 'MarkerSize', 8);
        
        % Σχεδίαση τελικού σημείου - ΜΑΥΡΟ Χ.
        final_x = x_hist(:, end);
        plot(final_x(1), final_x(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
        title({strategy_titles{j}, sprintf('Iters: %d, f=%.4f', k, f_hist(end))}, 'FontSize', 10);
        xlabel('x'); ylabel('y');
        hold off;

        fprintf('Method: %-15s | Iter: %3d | Final: (%.4f, %.4f) | f = %.4e\n', ...
            step_type, k, final_x(1), final_x(2), f_hist(end));
    end
    sgtitle(['Πορείες Σύγκλισης - Εκκίνηση: ' point_str]);
    
    % Plot: σύγκλιση της f.
    figure('Name', ['Σύγκλιση f(x) - Start: ' point_str], 'NumberTitle', 'off', 'Position', [100, 100, 1200, 350]);
    
    for j = 1:3
         step_type = step_strategies{j};
         [~, f_hist, ~] = steepest_descent(f, grad_f, x0, step_type, epsilon, max_iter);
         subplot(1, 3, j);
         plot(0:length(f_hist)-1, f_hist, '.-', 'Color', colors(j), 'LineWidth', 1.5, 'MarkerSize', 10);
         title(strategy_titles{j});
         xlabel('Επαναλήψεις (k)'); 
         ylabel('f(x_k)');
         grid on;
    end
    sgtitle(['Ρυθμός Σύγκλισης f(x) - Εκκίνηση: ' point_str]);
end