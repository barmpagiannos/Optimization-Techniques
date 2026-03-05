
% Τεχνικές Βελτιστοποίησης
% 1η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Αρχικό διάστημα (το δίνει η εκφώνηση).
a_init = -1;
b_init = 3;

% Τα ονόματα των τριών συναρτήσεων.
func_names = {'f_1(x)', 'f_2(x)', 'f_3(x)'};

% Οι παράγωγοι των συναρτήσεων που θα ελαχιστοποιήσουμε.
fp1 = @(x) log(5).*5.^x + (4 - 2*cos(x)).*sin(x);
fp2 = @(x) 2*(x-1)+exp(x-5).*(sin(x+3)+cos(x+3));
fp3 = @(x) -3*exp(-3*x)-2*cos(x-2).*(sin(x-2)-2);

% Πακετάρισμα των παραγώγων.
derivatives = {fp1, fp2, fp3};


%% Θέμα 4.1

l_vec_41 = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001];

eval_results_41 = zeros(length(derivatives), length(l_vec_41));

for i = 1:length(derivatives)
    fp_current = derivatives{i};
    for j = 1:length(l_vec_41)
        % Καλώ τη μέθοδο διχοτόμου με τη χρήση παραγώγου.
        [~, ~, ~, eval_count] = methodos_dixotomou_paragwgos(fp_current, a_init, b_init, l_vec_41(j));
        eval_results_41(i, j) = eval_count;
    end
end

% Δημιουργώ το γράφημα.
figure;
hold on;
for i = 1:length(derivatives)
    semilogx(l_vec_41, eval_results_41(i, :), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
end
hold off;
title('Θέμα 4.1 Υπολογισμοί παραγώγου f''(x) ως συνάρτηση του l');
xlabel('Τελικό εύρος l');
ylabel('Συνολικοί υπολογισμοί f''(x)');
legend(func_names, 'Location', 'best');
grid on;


%% Θέμα 4.2

% Διάφορες τιμές του l.
l_plot_vec = [0.1, 0.01, 0.001];

for i = 1:length(derivatives) 
    figure;
    hold on;
    title(['Θέμα 4.2 Σύγκλιση διαστήματος για ' func_names{i} '']);
    colors = ['r', 'g', 'b'];
    legend_entries = {};
    for j = 1:length(l_plot_vec)
        l_current = l_plot_vec(j);
        [~, ~, ~, ~, history] = methodos_dixotomou_paragwgos(derivatives{i}, a_init, b_init, l_current);
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
