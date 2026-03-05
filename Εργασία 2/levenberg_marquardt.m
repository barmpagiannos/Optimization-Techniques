
function [x_hist, f_hist, k] = levenberg_marquardt(f, grad_f, hess_f, x0, step_type, epsilon, max_iter)
    
    x_curr = x0;
    x_hist = x_curr;
    f_hist = f(x_curr);
    
    for k = 1:max_iter
        
        % Υπολογισμός κλίσης.
        g_k = grad_f(x_curr);
        
        % Έλεγχος τερματισμού.
        if norm(g_k) < epsilon
            break;
        end
        
        % Υπολογισμός Εσσιανού πίνακα.
        H_k = hess_f(x_curr);
        
        % Υπολογισμός του μ_k.
        % Βρίσκουμε τις ιδιοτιμές του Εσσιανού.
        e = eig(H_k);
        min_eigenval = min(e);
        
        % Αν η μικρότερη ιδιοτιμή είναι θετική, ο πίνακας είναι ήδη
        % θετικά ορισμένος. Δεν χρειάζεται διόρθωση (mu = 0 -> Newton).
        if min_eigenval > 0
            mu_k = 0;
        else
            % Αν είναι αρνητική ή μηδέν, θέτουμε mu_k λίγο μεγαλύτερο 
            % από την απόλυτη τιμή της, ώστε να γίνει θετική.
            % Προσθέτουμε ένα μικρό 'delta' (π.χ. 0.1) για ασφάλεια.
            mu_k = abs(min_eigenval) + 0.1; 
        end
        
        % Τροποποίηση Εσσιανού: H_mod = H + mu*I
        % Η εντολή eye(2) δημιουργεί τον μοναδιαίο πίνακα 2x2.
        H_mod = H_k + mu_k * eye(2);
        
        % Υπολογισμός κατεύθυνσης.
        d_k = -H_mod \ g_k;
        
        % Επιλογή βήματος.
        gamma_k = compute_step(f, grad_f, x_curr, d_k, step_type);
        
        % Ενημέρωση.
        x_next = x_curr + gamma_k * d_k;
        
        x_curr = x_next;
        x_hist = [x_hist, x_curr];
        f_hist = [f_hist, f(x_curr)];
    end
end
