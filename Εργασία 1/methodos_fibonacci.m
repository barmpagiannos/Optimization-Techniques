
function [a_final, b_final, k, eval_count, history] = methodos_fibonacci(f, a, b, l, epsilon)

L0 = b - a;
R = L0 / l;
F = [1, 1]; % F0 = F1 = 1.

while F(end) <= R
    F(end+1) = F(end) + F(end-1);
end
n = length(F) - 1; 
history.a(1) = a;
history.b(1) = b;

x1 = a + (F(n-1) / F(n+1)) * (b - a);
x2 = a + (F(n) / F(n+1)) * (b - a);

f1 = f(x1);
f2 = f(x2);
eval_count = 2;

for k = 1:(n-2)
    if f1 > f2
        a = x1;

        x1 = x2;
        f1 = f2;

        x2 = a + (F(n-k) / F(n-k+1)) * (b - a);

        if k < n-2
            f2 = f(x2);
            eval_count = eval_count + 1;
        end
        
    else
        b = x2;
   
        x2 = x1;
        f2 = f1;

        x1 = a + (F(n-k-1) / F(n-k+1)) * (b - a);
        if k < n-2
            f1 = f(x1);
            eval_count = eval_count + 1;
        end
    end
    history.a(k+1) = a;
    history.b(k+1) = b;
end 

x1_n = x1;
x2_n = x1_n + epsilon;

f1_n = f1;
f2_n = f(x2_n);
eval_count = eval_count + 1;

if f1_n > f2_n
    a = x1_n;
else
    b = x2_n; 
end
history.a(k+2) = a;
history.b(k+2) = b;
a_final = a;
b_final = b;
k = n - 1; 
end