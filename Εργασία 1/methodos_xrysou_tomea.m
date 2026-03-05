
function [a_final, b_final, k, eval_count, history] = methodos_xrysou_tomea(f, a, b, l)

gamma = (sqrt(5) - 1) / 2; % Σταθερά χρυσής αναλογίας = 0.61803.
k = 1;
% Αρχικό διάστημα [a1, b1]. Το αποθηκεύω.
history.a(k) = a;
history.b(k) = b;
x1 = a + (1 - gamma) * (b - a);
x2 = a + gamma * (b - a);
% Υπολογίζω τις τιμές f(x11) και f(x21).
f1 = f(x1);
f2 = f(x2);
eval_count = 2; 

while (b-a)>=l
    if f1 > f2 
        a = x1;
        x1 = x2; % το παλιό x2k γίνεται το νέο x1(k+1).
        f1 = f2; % η παλιά τιμή f2 γίνεται η νέα f1.
        x2 = a + gamma * (b - a);
        f2 = f(x2);
    else 
        b = x2;      
        x2 = x1; % το παλιό x1k γίνεται το νέο x2(k+1).
        f2 = f1; % η παλιά τιμή f1 γίνεται η νέα f2.
        x1 = a + (1 - gamma) * (b - a);
        f1 = f(x1);
    end
    eval_count = eval_count + 1;
    k = k + 1;
    history.a(k) = a;
    history.b(k) = b;
end
a_final = a;
b_final = b;
k = k - 1; % πόσες επαναλήψεις έγιναν.

end
