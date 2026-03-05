
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

functions = {f1, f2, f3};
func_names = {'f_1(x)', 'f_2(x)', 'f_3(x)'};

%% Θέμα 2.1

l_vec_21 = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001];

eval_results_21 = zeros(length(functions), length(l_vec_21));

% Καλώ τη μέθοδο χρυσού τομέα.
for i = 1:length(functions)
    f_current = functions{i};
    for j = 1:length(l_vec_21)
        [~, ~, ~, eval_count] = methodos_xrysou_tomea(f_current, a_init, b_init, l_vec_21(j));
        eval_results_21(i, j) = eval_count;
    end
end

% Τα βάζω όλα σε ένα γράφημα.
figure;
hold on;
for i = 1:length(functions)
    semilogx(l_vec_21, eval_results_21(i, :), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
end
hold off;
title('Θέμα 2.1  Υπολογισμοί της συνάρτησης f ως συνάρτηση του l');
xlabel('Τελικό εύρος l');
ylabel('Συνολικοί υπολογισμοί f(x)');
legend(func_names);
grid on;

%% Θέμα 2.2

l_plot_vec = [0.1, 0.01, 0.001];

for i = 1:length(functions) 
    figure; 
    hold on;
    title(['Θέμα 2.2 Σύγκλιση διαστήματος για ' func_names{i} '']);
    
    colors = ['r', 'g', 'b'];
    legend_entries = {};
    
    for j = 1:length(l_plot_vec)
        l_current = l_plot_vec(j);
        [~, ~, ~, ~, history] = methodos_xrysou_tomea(functions{i}, a_init, b_init, l_current);
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
