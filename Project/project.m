
% Τεχνικές Βελτιστοποίησης
% Project
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Δημιουργία Δεδομένων.
% Συνάρτηση παραγωγής δεδομένων (μόνο για εκπαίδευση και αξιολόγηση, όχι για χρήση στον GA).
true_f = @(u1, u2) sin(u1 + u2) .* sin(u2.^2);

num_train = 300;
% Πεδίο ορισμού για τις μεταβλητές u1, u2.
u1_min = -1; u1_max = 2;
u2_min = -2; u2_max = 1;

% Παράγουμε τα δεδομένα εκπαίδευσης.
u1_train = (u1_max - u1_min) * rand(num_train, 1) + u1_min;
u2_train = (u2_max - u2_min) * rand(num_train, 1) + u2_min;
y_train = true_f(u1_train, u2_train);

% Ρυθμίσεις Γενετικού Αλγορίθμου.
num_gaussians = 15; % Πλήθος Γκαουσιανών. Η εκφώνηση λέει το πολύ 15.
genes_per_gaussian = 5; % 5 παράμετροι της Γκαουσιανής: w, c1, σ1, c2, σ2.
chrom_len = num_gaussians * genes_per_gaussian; % Μήκος χρωμοσώματος.

pop_size = 60; % Πληθυσμός.
max_gen = 400; % Γενιές (δηλαδή επαναλήψεις).
mutation_rate = 0.05; % Ποσοστό μετάλλαξης (5%). Η πιθανότητα να αλλάξει τυχαία ένας αριθμός στο χρωμόσωμα.
mutation_step = 0.5; % Πόσο πολύ θα αλλάξει ο αριθμός αν συμβεί μετάλλαξη.

% Αρχικοποίηση Πληθυσμού (τυχαίες τιμές).
% c: [-2, 2], w: [-5, 5], sigma: [0.1, 2].
population = rand(pop_size, chrom_len);
% Κλιμάκωση (scaling): από [0, 1] -> [0, 4] -> [-2, 2].
population = (population * 4) - 2; 

% Αποθηκεύουμε το ιστορικό.
best_fitness_history = zeros(max_gen, 1);
best_solution = [];
min_error = inf;

fprintf('Χρησιμοποιώ %d Γκαουσιανές...\n', num_gaussians);

for gen = 1:max_gen
    
    % Αξιολόγηση.
    fitness = zeros(pop_size, 1);
    for i = 1:pop_size
        % Υπολογισμός σφάλματος (MSE) για κάθε χρωμόσωμα.
        y_pred = my_model(u1_train, u2_train, population(i, :), num_gaussians);
        mse = mean((y_train - y_pred).^2); % mean square error.
        fitness(i) = mse;
    end
    
    % Ο καλύτερος της γενιάς.
    [current_best_error, idx] = min(fitness);
    
    % Αποθήκευση ιστορικού.
    best_fitness_history(gen) = current_best_error;
    
    % Κρατάμε την καλύτερη λύση.
    if current_best_error < min_error
        min_error = current_best_error;
        best_solution = population(idx, :);
    end
    
    % Εκτύπωση ανά 10 γενιές.
    if mod(gen, 10) == 0 || gen == 1
        fprintf('Γενιά %d | Best MSE: %.5f\n', gen, current_best_error);
    end
    
    % Επιλογή.
    % Διαλέγουμε γονείς για την επόμενη γενιά.
    new_population = zeros(size(population));
    
    % Περνάμε τον καλύτερο στην επόμενη γενιά.
    new_population(1, :) = population(idx, :);
    start_idx = 2; % Οι υπόλοιποι θα παραχθούν.
    
    while start_idx <= pop_size
        % Διαλέγω 2 τυχαίους, κρατάω τον καλύτερο (Γονιός 1).
        r1 = randi(pop_size); r2 = randi(pop_size);
        if fitness(r1) < fitness(r2), p1 = population(r1,:); else, p1 = population(r2,:); end
        
        % Για Γονιό 2.
        r1 = randi(pop_size); r2 = randi(pop_size);
        if fitness(r1) < fitness(r2), p2 = population(r1,:); else, p2 = population(r2,:); end
        
        % Διασταύρωση (Crossover).
        % Arithmetic Crossover (Μέσος όρος με βάρος).
        alpha = rand();
        child = alpha * p1 + (1-alpha) * p2;
        
        % Μετάλλαξη (Mutation).
        % Για κάθε γονίδιο του παιδιού, υπάρχει πιθανότητα να αλλάξει λίγο.
        for g = 1:chrom_len
            if rand() < mutation_rate
                % Προσθέτουμε λίγο θόρυβο.
                child(g) = child(g) + mutation_step * randn(); 
            end
        end
        
        % Προσθήκη στον νέο πληθυσμό.
        new_population(start_idx, :) = child;
        start_idx = start_idx + 1;
    end
    
    population = new_population;
end

% Αποτέλεσμα.
fprintf('\nΤέλος Αλγορίθμου.\nΤελικό MSE: %.5e\n', min_error);

% Αξιολόγηση (Validation) σε ΝΕΑ δεδομένα, διαφορετικά από τα δεδομένα εκπαίδευσης.α
% Η εκφώνηση ζητάει να ελέγξουμε το μοντέλο σε διαφορετικό σύνολο δεδομένων.
fprintf('\n--- Έλεγχος Αξιολόγησης (Validation) ---\n');

% Δημιουργία 100 νέων τυχαίων σημείων (που δεν έχει ξαναδεί ο αλγόριθμος).
num_val = 100;
u1_val = (u1_max - u1_min) * rand(num_val, 1) + u1_min;
u2_val = (u2_max - u2_min) * rand(num_val, 1) + u2_min;
y_val_true = true_f(u1_val, u2_val); % Η πραγματική τιμή.

% Πρόβλεψη με το μοντέλο που βρήκαμε (best_solution).
y_val_pred = my_model(u1_val, u2_val, best_solution, num_gaussians);

% Υπολογισμός σφάλματος στα νέα δεδομένα.
mse_val = mean((y_val_true - y_val_pred).^2);

fprintf('Training MSE (Εκπαίδευση): %.5e\n', min_error);
fprintf('Validation MSE (Αξιολόγηση): %.5e\n', mse_val);

if mse_val < 0.01
    fprintf('Συμπέρασμα: Το μοντέλο γενικεύει ικανοποιητικά!\n');
else
    fprintf('Συμπέρασμα: Ίσως χρειάζονται περισσότερες γενιές ή Γκαουσιανές.\n');
end

% Γράφημα Σύγκλισης.
figure;
plot(best_fitness_history, 'LineWidth', 2);
title('Σύγκλιση MSE Γενετικού Αλγορίθμου');
xlabel('Γενιές');
ylabel('Mean Squared Error (MSE)');
grid on;

% 3D γράφημα.
figure;
u1_range = linspace(u1_min, u1_max, 40);
u2_range = linspace(u2_min, u2_max, 40);
[U1, U2] = meshgrid(u1_range, u2_range);

% Πραγματικό.
Y_true = true_f(U1, U2);
subplot(1,2,1);
surf(U1, U2, Y_true);
title('Real Function');
shading interp;

% Εκτίμηση (με το best_solution).
Y_pred = my_model(U1(:), U2(:), best_solution, num_gaussians);
Y_pred = reshape(Y_pred, size(U1));
subplot(1,2,2);
surf(U1, U2, Y_pred);
title('Estimated Model (GA)');
shading interp;

% Οι συναρτήσεις.
function y = my_model(u1, u2, chrom, K)
    y = zeros(size(u1));
    
    for i = 1:K
        idx = (i-1)*5;
        w  = chrom(idx+1);
        c1 = chrom(idx+2);
        s1 = chrom(idx+3);
        c2 = chrom(idx+4);
        s2 = chrom(idx+5);
        
        % Ασφάλεια: Το s (εύρος) δεν πρέπει να είναι 0 ή πολύ κοντά στο 0
        % ώστε να αποφευχθεί διαίρεση με το μηδέν.
        if abs(s1) < 1e-4, s1 = 1e-4; end
        if abs(s2) < 1e-4, s2 = 1e-4; end
        
        term = w .* exp( - ( (u1 - c1).^2 ./ (2*s1^2) + (u2 - c2).^2 ./ (2*s2^2) ) );
        y = y + term;
    end
end

