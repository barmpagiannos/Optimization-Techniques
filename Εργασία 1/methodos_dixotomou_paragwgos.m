
function [a_final, b_final, k, eval_count, history] = methodos_dixotomou_paragwgos(fp, a, b, l)

k = 1;
eval_count = 0; % Μετρητής υπολογισμών της f'(x).

history.a(k) = a;
history.b(k) = b;

while (b - a) >= l
    
    % Υπολογισμός του μέσου xk.
    xk = (a + b) / 2;
    
    % Υπολογισμός της παραγώγου στο xk.
    fp_val = fp(xk);
    eval_count = eval_count + 1;

    if abs(fp_val) < 1e-10
        % Το ελάχιστο βρέθηκε ακριβώς στο xk.
        a = xk;
        b = xk;
        break;
    elseif fp_val > 0
        % f'(x) > 0. Το ελάχιστο είναι αριστερά.
        % Νέο διάστημα αναζήτησης [a, xk].
        b = xk;
    else
        % f'(x) < 0. Το ελάχιστο είναι δεξιά.
        % Νέο διάστημα αναζήτησης (ΝΔΑ) [xk, b].
        a = xk;
    end
    k = k + 1; % k = k+1
    history.a(k) = a;
    history.b(k) = b;

end

a_final = a;
b_final = b;
k = k - 1;

end