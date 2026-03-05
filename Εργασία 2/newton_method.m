
function [x_hist, f_hist, k] = newton_method(f, grad_f, hess_f, x0, step_type, epsilon, max_iter)

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
        % Υπολογισμός εσσιανού πίνακα.
        H_k = hess_f(x_curr);
        % Υπολογισμός κατεύθυνσης Newton.
        % Λύνουμε το σύστημα H_k * d_k = -g_k.
        % Έλεγχος αν ο πίνακας είναι αντιστρέψιμος.
        if rcond(H_k) < 1e-12
            warning('Ο Εσσιανός πίνακας είναι σχεδόν ιδιάζων (singular) στο k=%d. Διακοπή.', k);
            break; 
        end
        
        d_k = -H_k \ g_k; 
        
        % Επιλογή βήματος.
        gamma_k = compute_step(f, grad_f, x_curr, d_k, step_type);
        
        % Ενημέρωση.
        x_next = x_curr + gamma_k * d_k;
        x_curr = x_next;
        x_hist = [x_hist, x_curr];
        f_hist = [f_hist, f(x_curr)];
    end
end
