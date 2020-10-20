function PureAbsorberSlab2
%1cm slab of pure absorber is discretized into n equal pieces. Discrete
%o source n the left side. Discrete angular discretization. We plot the 
%exiting angular flux

n = [25 50 100 200];
h = 1 ./ n;
L2_error_dd = zeros(1, length(n));
L2_error_sc = zeros(1, length(n));

for i = 1: length(n)
    %DEFINE SLAB AND ITS DISCRETIZATION
    edges = linspace(0, 1, n(i) + 1);
    slab = cell(1, n(i));
    for j = 1: n(i)
    slab{j} = material('pure absorber', edges(j), edges(j+1)); 
    end
    S_t = slab{1}.sig_t; %TOTAL CROSS-SECTION OF SLAB
    
    %DEFINE BOUNDARY CONDITIONS 
    psil = [0; 1];
    psir = [0; 0];
    %TYPE OF ANGULAR DISCRETIZATION
    flag = 'discrete';
    %TOLERANCE FOR CONVERGENCE
    tol = eps;
    
    %DIAMOND DIFFERENCE SOLUTION
    [~, ~, psil_dd] = diamond_difference(edges, slab, psil, psir, flag,...
    2, tol);
    %STEP CHARACTERISTICS SOLUTION
    [~, ~, psil_sc] = step_characteristics(edges, slab, psil, psir, flag,...
    2, tol);

    size(psil_dd)
    

    x = edges(2: n(i) + 1); %EXIT BOUNDARY 
    %TRUE SOLUTION
    sol = exp( - S_t * x); 
    size(sol)
    L2_error_dd(i) = sqrt(sum((sol - psil_dd(2, 2:(n(i) + 1))).^2)/n(i));
    L2_error_sc(i) = sqrt(sum((sol - psil_sc(2, 2:(n(i)+1))).^2 )/n(i));
end

figure
loglog(h, L2_error_dd, 'b.-',  h, L2_error_sc, 'r.-',...
    'MarkerSize', 20, 'LineWidth', 2)
legend({'Diamond difference', 'Step Characteristics'}, 'Location',...
 'west', 'FontSize', 14)
xlabel('mesh size')
ylabel('L^2 error')
title('Convergence rate in Pure Absorber of exiting angular flux',...
 'FontSize', 14)

end