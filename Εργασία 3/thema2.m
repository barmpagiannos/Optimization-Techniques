
% Τεχνικές Βελτιστοποίησης
% 3η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685
clc; clear; close all;

% Συνάρτηση f(x) = (1/3)x1^2 + 3x2^2
f = @(x) (1/3)*x(1)^2 + 3*x(2)^2;

% Κλίση (Gradient).
grad_f = @(x) [(2/3)*x(1); 6*x(2)];

% Περιορισμοί:
% -10 <= x1 <= 5,
% -8 <= x2 <= 12
x_min = [-10; -8];
x_max = [ 5; 12];

% Παράμετροι.
x0 = [5; -5]; % Σημείο εκκίνησης.
s_k = 5;
gamma_k = 0.1;
epsilon = 0.01;
max_iter = 50;

% Κλήση μεθόδου μέγιστης καθόδου με προβολή.
[x_hist, f_hist, k] = steepest_descent_projection(f, grad_f, x0, s_k, gamma_k, epsilon, max_iter, x_min, x_max);

% Πορεία σύγκλισης.
figure('Name', 'Thema 2: Projected Steepest Descent', 'NumberTitle', 'off');

% Ισοϋψείς καμπύλες.
[X, Y] = meshgrid(linspace(-12, 8, 100), linspace(-10, 14, 100));
Z = (1/3).*X.^2 + 3.*Y.^2;
contourf(X, Y, Z, 20, 'HandleVisibility', 'off'); hold on; colorbar;

% Σχεδίαση περιορισμών (box constraints).
rectangle('Position', [-10, -8, 15, 20], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--'); 
% Position: [x_min, y_min, width, height] -> [-10, -8, 5-(-10), 12-(-8)]

% Αρχικό και τελικό σημείο.
plot(x0(1), x0(2), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', sprintf('Start (%g, %g)', x0(1), x0(2)));
plot(x_hist(1, end), x_hist(2, end), 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Final');
plot(x_hist(1, :), x_hist(2, :), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Trajectory');

title(sprintf('Μέθοδος Προβολής (s_k=%g, \\gamma_k=%g)', s_k, gamma_k));
xlabel('x_1'); ylabel('x_2');
legend('show', 'Location', 'northeast');
grid on; hold off;

% Σύγκλιση της f(x).
figure('Name', 'Σύγκλιση της f(x)', 'NumberTitle', 'off');
plot(0:length(f_hist)-1, f_hist, 'b.-', 'LineWidth', 1.5);
title('Σύγκλιση Τιμής f(x)');
xlabel('Iterations (k)'); ylabel('f(x_k)');
grid on;

fprintf('Αποτελέσματα:\n');
fprintf('Επαναλήψεις: %d\n', k);
fprintf('Τελικό σημείο: (%.4f, %.4f)\n', x_hist(1,end), x_hist(2,end));
fprintf('Τελική τιμή f(x): %.4e\n', f_hist(end));

% Επίδειξη της προβολής. Βλέπουμε πως τραβάει το σημείο εντός του κουτιού.

x_min = [-10; -8];
x_max = [ 5; 12];

x_curr = x0; % Αρχικό Σημείο.

% Βήμα 1: Υπολογισμός Κλίσης
g = grad_f(x_curr);

% Βήμα 2: Το σημείο που θα βγει εκτός κουτιού.
w_k = x_curr - s_k * g; 

% Βήμα 3: Προβολή (bar_x_k) - Επαναφορά στο κουτί.
bar_x_k = min(max(w_k, x_min), x_max);

% Βήμα 4: Τελικό Σημείο (x_next).
x_next = x_curr + gamma_k * (bar_x_k - x_curr);

figure('Name', 'Mechanism of Projection', 'NumberTitle', 'off', 'Position', [100, 100, 900, 700]);

[X, Y] = meshgrid(linspace(-20, 15, 100), linspace(-10, 150, 100));
Z = (1/3).*X.^2 + 3.*Y.^2;

contour(X, Y, Z, 30, 'LineWidth', 1.0); hold on; 
colormap cool;

rectangle('Position', [-10, -8, 15, 20], 'EdgeColor', 'r', 'LineWidth', 3, 'LineStyle', '-');
text(-9, 10, 'Feasible Region', 'Color', 'r', 'FontWeight', 'bold');


% 1. Αρχικό Σημείο
plot(x_curr(1), x_curr(2), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(x_curr(1)+1, x_curr(2), 'Start (x_k)', 'FontSize', 10, 'FontWeight', 'bold');

% 2. Η Κίνηση ΕΚΤΟΣ Ορίων (x_k -> w_k).
quiver(x_curr(1), x_curr(2), w_k(1)-x_curr(1), w_k(2)-x_curr(2), 0, ...
    'k--', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
plot(w_k(1), w_k(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(w_k(1)+1, w_k(2), 'Unbounded Step (w_k)', 'Color', 'r', 'FontSize', 10);

% 3. Η ΠΡΟΒΟΛΗ (w_k -> bar_x_k).
quiver(w_k(1), w_k(2), bar_x_k(1)-w_k(1), bar_x_k(2)-w_k(2), 0, ...
    'm-', 'LineWidth', 2, 'MaxHeadSize', 0.5);
plot(bar_x_k(1), bar_x_k(2), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
text(bar_x_k(1)-6, bar_x_k(2)+2, 'Projection (\bar{x}_k)', 'Color', 'm', 'FontSize', 10, 'FontWeight', 'bold');

% 4. Η Τελική Κίνηση (x_k -> x_next).
plot([x_curr(1), bar_x_k(1)], [x_curr(2), bar_x_k(2)], 'g-', 'LineWidth', 2);
plot(x_next(1), x_next(2), 'bx', 'MarkerSize', 12, 'LineWidth', 3);
text(x_next(1)+1, x_next(2)-2, 'Next Point (x_{k+1})', 'Color', 'b', 'FontSize', 10);

title('Λεπτομέρεια 1ης Επανάληψης: Ο Μηχανισμός της Προβολής');
xlabel('x_1'); ylabel('x_2');
grid on; hold off;

% Κάνω focus εκεί που με ενδιαφέρει.
axis([-20 10 -15 160]);
