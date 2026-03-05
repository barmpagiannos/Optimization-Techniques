function yhat = model_predict(ind, U, M)
% U: Nx2, columns [u1 u2]
u1 = U(:,1);
u2 = U(:,2);

% Build features: each column is G_i(u1,u2)
N = size(U,1);
Phi = zeros(N, M);

for i = 1:M
    du1 = (u1 - ind.c1(i))./ind.s1(i);
    du2 = (u2 - ind.c2(i))./ind.s2(i);
    Phi(:,i) = exp(-0.5*(du1.^2 + du2.^2));
end

yhat = Phi * ind.w;
end
