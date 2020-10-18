function PureAbsorberSlab
%1cm slab of pure absorber discretized into n equal pieces. Discrete source
%on the left side

n = [10 100 1000 10000];
L2_error_dd = zeros(1, length(n));
L2_error_sc = zeros(1, length(n));

figdd1 = figure;
figdd2 = figure;
figsc1 = figure;
figsc2 = figure;
fig_conv = figure;

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
    [phi0_dd, ~] = diamond_difference(edges, slab, psil, psir, flag,...
    2, tol);
    %STEP CHARACTERISTICS SOLUTION
    [phi0_sc, ~] = step_characteristics(edges, slab, psil, psir, flag,...
    2, tol);
    
    x = edges(1: n(i)) + 1/(2 * n(i));%DISCRETIZED DOMAIN
    %TRUE SOLUTION
    sol = exp( - S_t * x); 
    
    L2_error_dd(i) = sqrt(sum((sol - phi0_dd).^2)/n(i));
    L2_error_sc(i) = sqrt(sum((sol - phi0_sc).^2 )/n(i));
    
    figure(figdd1)
    sgtitle('Diamond difference compared to True Solution in Pure Absorber')
    subplot(2, 2, i)
    plot(x, phi0_dd, 'r.', x, sol,'k', 'MarkerSize', 10, 'LineWidth', 2 )
    xlabel('x')
    ylabel('\phi_0')
    str = sprintf('mesh size = %.4f', 1/n(i));
    title(str)
    
    figure(figdd2)
    sgtitle('Diamond difference relative error in Pure Absorber')
    subplot(2, 2, i)
    plot(x, abs(phi0_dd - sol)./sol, 'b', 'LineWidth', 2 )
    xlabel('x')
    ylabel('$\mathbf{\frac{|\phi_0 - \psi_{exact}|}{\psi_{exact}}}$',...
        'Interpreter', 'latex', 'FontSize', 25)
    str = sprintf('mesh size = %.4f', 1/n(i));
    title(str)
      
    figure(figsc1)
    sgtitle('Step characteristics compared to True Solution in Pure Absorber')
    subplot(2, 2, i)
    plot(x, phi0_sc, 'r.', x, sol,'k', 'MarkerSize', 10, 'LineWidth', 2 )
    xlabel('x')
    ylabel('\phi_0')
    str = sprintf('mesh size = %.4f', 1/n(i));
    title(str)
    
    figure(figsc2)
    sgtitle('Step characteristics relative error in Pure Absorber')
    subplot(2, 2, i)
    plot(x, abs(phi0_sc - sol)./sol, 'r', 'LineWidth', 2 )
    xlabel('x')
    ylabel('$\mathbf{\frac{|\phi_0 - \psi_{exact}|}{\psi_{exact}}}$',...
        'Interpreter', 'latex', 'FontSize', 25)
    str = sprintf('mesh size = %.4f', 1/n(i));
    title(str)
end
figure(fig_conv)
loglog(1./n, L2_error_dd,'b.-',  1./n, L2_error_sc, 'r.-',...
    'MarkerSize', 20, 'LineWidth', 2)
legend('Diamond difference', 'Step Characteristic')
xlabel('mesh size')
ylabel('L^2 error')
title('Convergence rate in Pure Absorber')



end