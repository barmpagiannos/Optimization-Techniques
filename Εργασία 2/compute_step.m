
% Με αυτή τη βοηθητική συνάρτηση υπολογίζουμε το βήμα γ_k με βάση τους
% τρεις τρόπους που ζητάει η εργασία:
% α) σταθερό
% β) τέτοιο ώστε να ελαχιστοποιεί την f(x_k+γ_kd_k)
% γ) βάσει του κανόνα Armijio.

function gamma = compute_step(f, grad_f, x_k, d_k, type)

    switch type
        case 'constant'
            % α) Σταθερό βήμα (επιλέγουμε μια μικρή τιμή, της επιλογής μας).
            gamma = 0.1; 
            
        case 'min'
            % β) Ελαχιστοποίηση της f(x_k + gamma * d_k).
            phi = @(g) f(x_k + g * d_k);
            % Χρήση μικρότερου άνω ορίου για να αποφύγουμε "μεγάλα άλματα"
            % που βγάζουν έξω από την περιοχή ενδιαφέροντος.
            lower_bound = 0;
            upper_bound = 1;
            gamma = fminbnd(phi, lower_bound, upper_bound);
            % Για ασφάλεια: αν επιστραφεί μεγάλο gamma το περιορίζουμε.
            gamma = min(gamma, 2);
            
        case 'armijo'
            % γ) Κανόνας Armijo.
            alpha = 1e-4;
            beta = 0.5;
            s = 1; % Αρχικό βήμα δοκιμής.
            gamma = s;
            m_k = 0;
            g_k = grad_f(x_k);
            % Σιγουρευόμαστε ότι d_k είναι κατεύθυνση καθόδου:
            if g_k' * d_k >= 0
                % Αν δεν είναι καθόδου, αντιστρέφουμε το πρόσημο του d_k.
                d_k = -d_k;
            end
            while f(x_k + gamma * d_k) > f(x_k) + alpha * gamma * (g_k' * d_k)
                m_k = m_k + 1;
                gamma = s * (beta^m_k);
                if m_k > 40
                    break; 
                end
            end
    end
end
