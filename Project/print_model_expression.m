function print_model_expression(ind, M, tau_w)
% Prints a compact analytic expression by keeping only active Gaussians
idx = find(abs(ind.w) > tau_w);
fprintf('\n--- Estimated analytic expression (compact) ---\n');
fprintf('f_hat(u1,u2) = sum_{i in active} w_i * exp(-( (u1-c1_i)^2/(2*s1_i^2) + (u2-c2_i)^2/(2*s2_i^2) ))\n');
fprintf('Active terms (|w|>%.3f): %d\n\n', tau_w, numel(idx));

for k = 1:numel(idx)
    i = idx(k);
    fprintf(['Term %2d:  w= %+ .6f,  c1= %+ .6f,  c2= %+ .6f,  s1= %.6f,  s2= %.6f\n'], ...
        i, ind.w(i), ind.c1(i), ind.c2(i), ind.s1(i), ind.s2(i));
end
fprintf('---------------------------------------------\n');
end
