
function [a_final, b_final, k, eval_count, history] = methodos_dixotomou(f, a, b, l, epsilon)

k = 1; % αρχική τιμή του k
eval_count = 0;
history.a(k) = a; % κρατάμε ιστορικό με τις τιμές a και b.
history.b(k) = b;
while (b-a)>=l
    x1k = (a+b)/2-epsilon;
    x2k = (a+b)/2+epsilon;
    f1 = f(x1k);
    f2 = f(x2k);
    eval_count = eval_count+2;
    if f1<f2
        b = x2k; 
    else
        a = x1k;
    end
    % Αύξηση μετρητή επανάληψης.
    k = k + 1; % Θέτουμε k = k+1.
    history.a(k) = a;
    history.b(k) = b;
end
% Τέλος αλγορίθμου. Το ελάχιστο ανήκει στο [a, b].
a_final = a;
b_final = b;
k = k - 1; 

end
