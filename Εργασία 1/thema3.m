
% Τεχνικές Βελτιστοποίησης
% 1η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Αρχικό διάστημα (το δίνει η εκφώνηση).
a_init = -1;
b_init = 3;

% Η σταθερά 'epsilon' που χρειάζεται η Fibonacci στο Βήμα 5 (έτσι όπως το λέει το βιβλίο).
% Επιλέγουμε μια πολύ μικρή τιμή, δεν επηρεάζει τον αριθμό βημάτων n.
epsilon = 1e-6; 

% Οι συναρτήσεις που θα ελαχιστοποιηθούν.
f1 = @(x) 5^x+(2-cos(x)).^2;
f2 = @(x) (x-1).^2+exp(x-5).*sin(x+3);
f3 = @(x) exp(-3*x)-(sin(x-2)-2).^2;

functions = {f1, f2, f3};
func_names = {'f_1(x)', 'f_2(x)', 'f_3(x)'};


%% Θέμα 3.1

l_vec_31 = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001];

eval_results_31 = zeros(length(functions), length(l_vec_31));

% Καλώ τη μέθοδο Fibonacci.
for i = 1:length(functions)
    f_current = functions{i};
    for j = 1:length(l_vec_31)
        [~, ~, ~, eval_count] = methodos_fibonacci(f_current, a_init, b_init, l_vec_31(j), epsilon);
        eval_results_31(i, j) = eval_count;
    end
end

% Τα βάζω όλα σε ένα γράφημα.
figure;
hold on;
for i = 1:length(functions)
    semilogx(l_vec_31, eval_results_31(i, :), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
end
hold off;
title('Θέμα 3.1 Υπολογισμοί της συνάρτησης f ως συνάρτηση του l');
xlabel('Τελικό Εύρος l');
ylabel('Συνολικοί Υπολογισμοί f(x)');
legend(func_names, 'Location', 'best');
grid on;


%% Θέμα 3.2

% Διάφορες τιμές του l.
l_plot_vec = [0.1, 0.01, 0.001];

% 3 διαγράμματα, ένα για κάθε συνάρτηση.
for i = 1:length(functions)
    figure;
    hold on;
    title(['Θέμα 3.2 Σύγκλιση διαστήματος για ' func_names{i} '']);
    colors = ['r', 'g', 'b'];
    legend_entries = {};
    
    for j = 1:length(l_plot_vec)
        l_current = l_plot_vec(j);
        % Καλώ τη μέθοδο Fibonacci.
        [~, ~, ~, ~, history] = methodos_fibonacci(functions{i}, a_init, b_init, l_current, epsilon);
        k_axis = 0:length(history.a)-1;
        plot(k_axis, history.a, '--', 'Color', colors(j), 'LineWidth', 1.5);
        plot(k_axis, history.b, '-', 'Color', colors(j), 'LineWidth', 1.5);
        legend_entries{end+1} = ['a_k (l=' num2str(l_current) ')'];
        legend_entries{end+1} = ['b_k (l=' num2str(l_current) ')'];
    end
    hold off;
    xlabel('Επανάληψη k');
    ylabel('Τιμή άκρων διαστήματος [a_k, b_k]');
    legend(legend_entries, 'Location', 'best');
    grid on;
end
