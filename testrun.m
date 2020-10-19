
n = 504;
edges = linspace(0, 1, n+1);
slab = cell(1, n);
for i = 1: n
    slab{i} = material('fuel', edges(i), edges(i+1)); 
end
flag = 'glq';
k = 12;
psi_l = 1 * ones(k,1);
psi_r = 1 * ones(k,1);

tol = 1e-06;

[phi0, phi1] = step_characteristics(edges, slab, psi_l, psi_r, flag,...
    k, tol);

%disp(iter)

x = edges(1:length(edges)-1) + (edges(2) - edges(1)) / 2;
phi00 = n * (1 - exp(-sqrt(0.6)/n))/ (sqrt(0.6)); %mean([1, exp(-x(2) * sqrt(0.6))]) ;

disp(phi0(1))
figure
plot(x, (phi0), 'k.');%, x, phi00 * exp(-x * sqrt(0.6)), 'r', 'LineWidth', 2)
figure
plot(x, (phi1), 'k.');%, x, .129099 * phi00 * exp(-x * sqrt(0.6)), 'r', 'LineWidth', 2)
%, edges(1:length(edges) - 1),...
    % -exp(- 10 * edges(1:length(edges)-1)), 'k')

%