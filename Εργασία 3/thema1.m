
% Τεχνικές Βελτιστοποίησης
% 3η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc; clear; close all;

% Συνάρτηση f(x) = (1/3)x1^2 + 3x2^2
f = @(x) (1/3)*x(1)^2 + 3*x(2)^2;

% Κλίση (Gradient).
grad_f = @(x) [
    (2/3)*x(1); % df/dx1.
    6*x(2) % df/dx2.
];

x0 = [5; -5]; % Αρχικό σημείο.
epsilon = 0.001; % ε = 0.001.
max_iter = 50; % μέγιστες επαναλήψεις.
gammas = [0.1, 0.3, 3, 5]; % τα 4 βήματα που ζητάει η άσκηση.

% Ορισμός πεδίου ορισμού.
% Δημιουργούμε ένα πλέγμα (grid) σημείων (x,y) για να σχεδιάσουμε τη συνάρτηση.
% Επιλέγουμε το διάστημα [-12, 12] που καλύπτει τα σημεία ενδιαφέροντος.
[X, Y] = meshgrid(linspace(-12, 12, 100), linspace(-12, 12, 100));
Z = (1/3).*X.^2 + 3.*Y.^2;

% 3D γράφημα της f(x).
figure('Name', 'Figure 1: 3D Surface Plot', 'NumberTitle', 'off');
surf(X, Y, Z);
shading interp;
colormap jet;
colorbar;
xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
title('3D Γράφημα της f(x) = 1/3 x_1^2 + 3 x_2^2');
view(45, 30);
grid on;

% Ισοϋψείς καμπύλες της f(x).
figure('Name', 'Figure 2: Contour Plot', 'NumberTitle', 'off');
contourf(X, Y, Z, 20, 'HandleVisibility', 'off'); 
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Ισοϋψείς Καμπύλες της f(x)');
grid on; hold on;

% Αρχικό σημείο και σημείο ελαχίστου.
plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Ελάχιστο (0,0)');
plot(5, -5, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', 'Αρχή (5,-5)');

legend('show', 'Location', 'northeast');
hold off;

% Figures για τις πορείες σύγκλισης και για τη σύγκλιση της f(x).
figure('Name', 'Figure 3: Πορείες Σύγκλισης', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);
figure('Name', 'Figure 4: Σύγκλιση Τιμής f(x)', 'NumberTitle', 'off', 'Position', [150, 150, 1000, 800]);

% Τρέχω τη μέθοδο μέγιστης καθόδου για κάθε τιμή του gamma.
for i = 1:4
    gamma_val = gammas(i);
    
    % Εκτέλεση της μεθόδου.
    [x_hist, f_hist, k] = steepest_descent(f, grad_f, x0, gamma_val, epsilon, max_iter);
    
    set(0, 'CurrentFigure', 3);
    subplot(2, 2, i);
    
    contourf(X, Y, Z, 20, 'HandleVisibility', 'off'); hold on; 
    
    plot(x0(1), x0(2), 'kp', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
    plot(x_hist(1, end), x_hist(2, end), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Πορεία σύγκλισης.
    plot(x_hist(1, :), x_hist(2, :), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8);
    
    title(sprintf('Βήμα \\gamma = %g (Iters: %d)', gamma_val, k));
    xlabel('x_1'); ylabel('x_2');
    grid on; hold off;
    
    set(0, 'CurrentFigure', 4);
    subplot(2, 2, i);
    
    plot(0:length(f_hist)-1, f_hist, 'b.-', 'LineWidth', 1.5);
    
    title(sprintf('Σύγκλιση f(x) με \\gamma = %g', gamma_val));
    xlabel('Επαναλήψεις (k)'); ylabel('f(x_k)');
    grid on;
    
    fprintf('Gamma: %.1f | Iter: %3d | Final f: %.4e | Point: (%.2f, %.2f)\n', ...
        gamma_val, k, f_hist(end), x_hist(1,end), x_hist(2,end));
end

%% Μέθοδος μέγιστης καθόδου.
function [x_hist, f_hist, k] = steepest_descent(f, grad_f, x0, gamma, epsilon, max_iter)
    
    x_curr = x0; % Τρέχον σημείο.
    x_hist = x_curr; % Αποθήκευση ιστορικού.
    f_hist = f(x_curr);
    
    for k = 1:max_iter
        % Βήμα 1. Υπολογισμός κλίσης:
        g_k = grad_f(x_curr);
        
        % Βήμα 2. Έλεγχος τερματισμού:
        if norm(g_k) < epsilon
            break;
        end
        
        % Βήμα 3. Κατεύθυνση αναζήτησης: d_k = -grad f(x_k).
        d_k = -g_k;
        
        % Βήμα 4. Σταθερό βήμα.
        x_next = x_curr + gamma * d_k;
        
        % Αποθήκευση για την κατασκευή των γραφημάτων μετά.
        x_curr = x_next;
        x_hist = [x_hist, x_curr]; % Προσθήκη στο ιστορικό.
        f_hist = [f_hist, f(x_curr)];
    end
end
