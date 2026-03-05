
% Τεχνικές Βελτιστοποίησης
% Project
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Δημιουργία Δεδομένων.
% Συνάρτηση παραγωγής δεδομένων (μόνο για εκπαίδευση και αξιολόγηση, όχι για χρήση στον GA).
true_f = @(u1, u2) sin(u1 + u2) .* sin(u2.^2);

% Δεδομένα εκπαίδευσης.
num_train = 300;
% Πεδίο ορισμού για τις μεταβλητές u1, u2.
u1_min = -1; u1_max = 2;
u2_min = -2; u2_max = 1;
% Παράγουμε τα δεδομένα εκπαίδευσης.
u1_train = (u1_max - u1_min) * rand(num_train, 1) + u1_min;
u2_train = (u2_max - u2_min) * rand(num_train, 1) + u2_min;
y_train = true_f(u1_train, u2_train);

% Δεδομένα αξιολόγησης (διαφορετικά από τα δεδομένα εκπαίδευσης).
num_val = 100;
u1_val = (u1_max - u1_min) * rand(num_val, 1) + u1_min;
u2_val = (u2_max - u2_min) * rand(num_val, 1) + u2_min;
y_val_true = true_f(u1_val, u2_val);

% Ρυθμίσεις γενετικού αλγορίθμου.
num_gaussians = 15;
genes_per_gaussian = 5; % w, c1, sigma1, c2, sigma2.
chrom_len = num_gaussians * genes_per_gaussian;

% Παράμετροι.
pop_size = 100;
max_gen = 500;
crossover_prob = 0.8;
mutation_rate = 0.01;

% Αρχικοποίηση Πληθυσμού.
population = (rand(pop_size, chrom_len) * 4) - 2;

best_fitness_history = zeros(max_gen, 1);
min_mse_history = zeros(max_gen, 1); % Για να βλέπουμε το σφάλμα.
best_solution = [];
global_min_mse = inf;

fprintf('Εκκίνηση ΓΑ (Τροχός της τύχης και διασταύρωση μονού σημείου)\n');

for gen = 1:max_gen
    
    % Υπολογισμός MSE και fitness.
    mse_values = zeros(pop_size, 1);
    for i = 1:pop_size
        y_pred = my_model(u1_train, u2_train, population(i, :), num_gaussians);
        mse_values(i) = mean((y_train - y_pred).^2);
    end
    
    % Κρατάμε το καλύτερο.
    [current_min_mse, best_idx] = min(mse_values);
    if current_min_mse < global_min_mse
        global_min_mse = current_min_mse;
        best_solution = population(best_idx, :);
    end
    min_mse_history(gen) = current_min_mse;
    
    % Μετατροπή MSE σε fitness.
    raw_fitness = 1 ./ (mse_values + 1e-6);
    
    % Κλιμακώνουμε την ικανότητα για το λόγο που είπαμε στην αναφορά.
    f_min = min(raw_fitness);
    f_max = max(raw_fitness);
    if f_max > f_min
        scaled_fitness = (raw_fitness - f_min) / (f_max - f_min);
    else
        scaled_fitness = ones(pop_size, 1);
    end
    
    total_fitness = sum(scaled_fitness);
    if total_fitness == 0, probs = ones(pop_size,1)/pop_size; else, probs = scaled_fitness / total_fitness; end
    cum_probs = cumsum(probs);
    
    new_population = zeros(size(population));
    
    % Κρατάμε τον καλύτερο.
    new_population(1, :) = population(best_idx, :);
    start_k = 2;
    
    % Παραγωγή νέας γενιάς ανά ζεύγη.
    while start_k <= pop_size
        % Επιλογή Γονέα 1.
        r = rand();
        parent1_idx = find(cum_probs >= r, 1, 'first');
        if isempty(parent1_idx), parent1_idx = pop_size; end
        
        % Επιλογή Γονέα 2.
        r = rand();
        parent2_idx = find(cum_probs >= r, 1, 'first');
        if isempty(parent2_idx), parent2_idx = pop_size; end
        
        p1 = population(parent1_idx, :);
        p2 = population(parent2_idx, :);
        
        % Πιθανότητα διασταύρωσης.
        if rand() < crossover_prob
            % Τυχαίο σημείο τομής (από 1 έως chrom_len-1).
            c_point = randi([1, chrom_len-1]);
            
            % Ανταλλαγή τμημάτων.
            child1 = [p1(1:c_point), p2(c_point+1:end)];
            child2 = [p2(1:c_point), p1(c_point+1:end)];
        else
            child1 = p1;
            child2 = p2;
        end
        
        % Μετάλλαξη.
        % Εφαρμόζεται με πολύ μικρή πιθανότητα σε κάθε γονίδιο
        % Εδώ προσθέτουμε θόρυβο αντί για αντιστροφή ψηφίου.
        for g = 1:chrom_len
            if rand() < mutation_rate
                child1(g) = child1(g) + randn(); 
            end
            if rand() < mutation_rate
                child2(g) = child2(g) + randn();
            end
        end
        
        % Αποθήκευση απογόνων.
        new_population(start_k, :) = child1;
        if start_k + 1 <= pop_size
            new_population(start_k + 1, :) = child2;
        end
        start_k = start_k + 2;
    end
    
    population = new_population;
    
    % Εκτύπωση ανά 50 γενιές.
    if mod(gen, 50) == 0 || gen == 1
        fprintf('Γενιά %d | Best MSE: %.5e\n', gen, global_min_mse);
    end
end

fprintf('\n--- ΤΕΛΟΣ ΑΛΓΟΡΙΘΜΟΥ ---\n');
fprintf('Training MSE: %.5e\n', global_min_mse);

% Έλεγχος στο validation set.
y_val_pred = my_model(u1_val, u2_val, best_solution, num_gaussians);
val_mse = mean((y_val_true - y_val_pred).^2);
fprintf('Validation MSE: %.5e\n', val_mse);

figure('Name', 'Convergence');
semilogy(min_mse_history, 'LineWidth', 2);
title('Ιστορικό Σφάλματος (MSE)');
xlabel('Generations'); ylabel('MSE (log scale)');
grid on;

figure('Name', 'Model Comparison');
% Πλέγμα για 3D απεικόνιση.
[U1, U2] = meshgrid(linspace(-1, 2, 40), linspace(-2, 1, 40));
Y_real = true_f(U1, U2);
Y_est = reshape(my_model(U1(:), U2(:), best_solution, num_gaussians), size(U1));

subplot(1,2,1);
surf(U1, U2, Y_real);
title('Πραγματική (Target)');
shading interp;
axis tight;
subplot(1,2,2);
surf(U1, U2, Y_est);
title('Προσέγγιση (GA)');
shading interp;
axis tight;


%% Απαραίτητες συναρτήσεις.
function y = my_model(u1, u2, chrom, K)
    y = zeros(size(u1));
    for k = 1:K
        idx = (k-1)*5;
        w  = chrom(idx+1);
        c1 = chrom(idx+2);
        s1 = chrom(idx+3);
        c2 = chrom(idx+4);
        s2 = chrom(idx+5);
        
        % Αποφυγή διαίρεσης με το 0.
        if abs(s1) < 1e-6, s1 = 1e-6; end
        if abs(s2) < 1e-6, s2 = 1e-6; end
        
        term = w .* exp( - ( (u1 - c1).^2 ./ (2*s1^2) + (u2 - c2).^2 ./ (2*s2^2) ) );
        y = y + term;
    end
end
