
% Τεχνικές Βελτιστοποίησης
% 1η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Αρχικό διάστημα (το δίνει η εκφώνηση).
a_init = -1;
b_init = 3;

% Οι συναρτήσεις που θα ελαχιστοποιηθούν.
f1 = @(x) 5^x+(2-cos(x)).^2;
f2 = @(x) (x-1).^2+exp(x-5).*sin(x+3);
f3 = @(x) exp(-3*x)-(sin(x-2)-2).^2;

% Τις πακετάρω και βάζω και ονόματα. Θα βοηθήσουν μετά στα plots.
functions = {f1, f2, f3};
func_names = {'f_1(x)', 'f_2(x)', 'f_3(x)'};

%% Θέμα 1.1

l_fixed = 0.01; % Κρατάω το l σταθερό.

% Χρησιμοποιούμε ένα διάνυσμα τιμών epsilon.
% Σιγουρευόμαστε ότι όλες οι τιμές είναι < l/2 = 0.005, διότι πρέπει 
% l >= 2e ώστε να τερματίσει ο αλγόριθμος.
epsilon_vec_11 = [1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3, 4e-3]; % 11 = θέμα 1, ερώτημα 1.

% Σε αυτό το πίνακα θα βάλουμε τα αποτελέσματα.
eval_results_11 = zeros(length(functions), length(epsilon_vec_11));

% Καλώ τη μέθοδο διχοτόμησης.
for i = 1:length(functions)
    f_current = functions{i};
    for j = 1:length(epsilon_vec_11)
        [~, ~, ~, eval_count] = methodos_dixotomou(f_current, a_init, b_init, l_fixed, epsilon_vec_11(j));
        eval_results_11(i, j) = eval_count;
    end
end

% Τα βάζω όλα σε ένα γράφημα.
figure;
hold on;
for i = 1:length(functions)
    semilogx(epsilon_vec_11, eval_results_11(i, :), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
end
hold off;
title('Θέμα 1.1 Υπολογισμοί συνάρτησης f ως συνάρτηση του ε');
xlabel('Τιμή ε (λογαριθμική κλίμακα)');
ylabel('Συνολικοί υπολογισμοί f(x)');
legend(func_names);
grid on;

%% Θέμα 1.2

epsilon_fixed = 0.001; % Κρατάω το epsilon σταθερό.

l_vec_12 = [0.003, 0.005, 0.01, 0.05, 0.1, 0.5]; % 12: θέμα 1, ερώτημα 2.

eval_results_12 = zeros(length(functions), length(l_vec_12));

% Καλώ τη μέθοδο διχοτόμησης.
for i = 1:length(functions) 
    f_current = functions{i};
    for j = 1:length(l_vec_12)
        [~, ~, ~, eval_count] = methodos_dixotomou(f_current, a_init, b_init, l_vec_12(j), epsilon_fixed);
        eval_results_12(i, j) = eval_count;
    end
end

% Τα βάζω όλα σε ένα plot.
figure;
hold on;
for i = 1:length(functions)
    semilogx(l_vec_12, eval_results_12(i, :), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
end
hold off;
title('Θέμα 1.2 Υπολογισμοί συνάρτησης f ως συνάρτηση του l');
xlabel('Τελικό εύρος l (λογαριθμική κλίμακα)');
ylabel('Συνολικοί υπολογισμοί f(x)');
legend(func_names);
grid on;

%% Θέμα 1.3

l_plot_vec = [0.1, 0.01, 0.001]; % διάφορες τιμές του τελικού εύρους αναζήτησης.
epsilon_plot = 0.0001; % διαλέγω μια τιμή για το epsilon.

colors = ['r', 'g', 'b'];

for i = 1:length(functions) 
    for j = 1:length(l_plot_vec) % Βρόχος για l=0.1, 0.01, 0.001.
        l_current = l_plot_vec(j);
        current_color = colors(j);
        figure;
        % Καλώ τη μέθοδο της διχοτόμου.
        [~, ~, ~, ~, history] = methodos_dixotomou(functions{i}, a_init, b_init, l_current, epsilon_plot);
        k_axis = 0:length(history.a)-1;
        plot(k_axis, history.a, '--', 'Color', current_color, 'LineWidth', 1.5);
        hold on;
        plot(k_axis, history.b, '-', 'Color', current_color, 'LineWidth', 1.5);
        hold off;
        title_text = sprintf('Θέμα 1.3 Σύγκλιση για %s (με l = %g)', func_names{i}, l_current);
        title(title_text);
        xlabel('Επανάληψη k');
        ylabel('Τιμή άκρων διαστήματος [a_k, b_k]');
        legend_entries = {['a_k (l=' num2str(l_current) ')'], ['b_k (l=' num2str(l_current) ')']};
        legend(legend_entries, 'Location', 'best');
        grid on;
    end
end
