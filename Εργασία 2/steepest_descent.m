
% Αυτή η συνάρτηση παίρνει το αρχικό σημείο και το βήμα και επιστρέφει την
% πορεία σύγκλισης. Υλοποιεί την μέθοδο της μέγιστης καθόδου.
function [x_hist, f_hist, k] = steepest_descent(f, grad_f, x0, step_type, epsilon, max_iter)

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
        
        % Βήμα 4. Επιλογή βήματος:
        % Κλήση της βοηθητικής συνάρτησης compute_step για τον υπολογισμό
        % του βήματος.
        gamma_k = compute_step(f, grad_f, x_curr, d_k, step_type);
        
        % Βήμα 5. Ενημέρωση σημείου:
        x_next = x_curr + gamma_k * d_k;
        
        % Αποθήκευση για την κατασκευή των γραφημάτων μετά.
        x_curr = x_next;
        x_hist = [x_hist, x_curr]; % Προσθήκη στο ιστορικό.
        f_hist = [f_hist, f(x_curr)];
    end
end