function PureScatteringSlab
%1cm slab of scatterer material discretized into n equal pieces. Isotropic
%source on left boundary. Gauss-Legendre quadrature for angular
%discretization n = 2,4, 6, 8

n = [25 50 100 200]; %SPATIAL DISCRETIZATIONS

S_t = material('scatterer').sig_t; %TOTAL CROSS-SECTION OF SLAB
S_s = material('scatterer').sig_s0; %SCATTERING CROSS-SECTION OF SLAB
S_a = S_t - S_s; %ABSORPTION CROSS-SECTION
D = 1 / (3 * S_t); %DIFFUSION CONSTANT
A = sqrt(S_a / D);  %will be used in diffusion solution

flag = 'glq'; %GAUSS LEGENDRE ANGULAR DISCRETIZATION
m = [2 4 6 8]; %NUMBER OF GAUSS LEGENDRE QUADRATURE POINTS

%BOUNDARY CONDITIONS
psil = cell(1, length(m));
psir = cell(1, length(m));
for j = 1 : length(m) 
    psil{j} = ones(m(j), 1); % isotropic source on left
    psir{j} = zeros(m(j), 1);% vacuum on right
end

%TOLERANCE FOR CONVERGENCE
tol = 1e-06; 

phi0diff = cell(length(n), 1); % 0th diffusion solution GOES HERE
phi1diff = cell(length(n), 1); % 1st diffusion solution GOES HERE
  
phi0dd = cell(length(n), length(m)); %for 0th diamond difference solution
phi1dd = cell(length(n), length(m)); %for 1st dd solution
phi0sc = cell(length(n), length(m)); %for 0th step characteristics solution
phi1sc = cell(length(n), length(m)); %for 1st sc solution
     
%L-2 error compared diffusion 
err0dd = cell(1, length(m));%cell(length(n), length(m));
err1dd = cell(1, length(m));%cell(length(n), length(m)); 
err0sc = cell(1, length(m));%cell(length(n), length(m));
err1sc = cell(1, length(m));% cell(length(n), length(m));

for k = 1: length(m)
    err0dd{k} = zeros(1, length(n));
    err1dd{k} = zeros(1, length(n));
    err0sc{k} = zeros(1, length(n));
    err1sc{k} = zeros(1, length(n));
end

for i = 1: length(n)
    %SLAB DISCRETIZATION
    edges = linspace(0, 10, n(i) + 1);
    h = abs(edges(2) - edges(1)); %mesh size
   
    %DEFINE THE SLAB
    slab = cell(1, n(i));
    for j = 1: n(i)
    slab{j} = material('scatterer', edges(j), edges(j+1)); 
    end
    
    z = edges(1: n(i)) +  h / 2  ;%DOMAIN POINTS
    
    %FILL IN DIFFUSION SOLUTIONS
     phi0diff{i} = 1.55 * exp(- z * A);
     phi1diff{i} = 0.2 * exp(- z * A);
     
     
    for k = 1 : length(m) 
        [phi0dd{i, k}, phi1dd{i, k}] = diamond_difference(edges, slab,...
            psil{k}, psir{k}, flag, m(k), tol) ;
        
        [phi0sc{i, k}, phi1sc{i, k}] = step_characteristics(edges, slab,...
            psil{k}, psir{k}, flag, m(k), tol) ;
        
        err0dd{k}(i) = sqrt(sum((phi0diff{i} - phi0dd{i,k}).^2)* h);
        err1dd{k}(i) = sqrt(sum((phi1diff{i} - phi1dd{i,k}).^2)* h);
        err0sc{k}(i) = sqrt(sum((phi0diff{i} - phi0sc{i,k}).^2)* h);
        err1sc{k}(i) = sqrt(sum((phi1diff{i} - phi1sc{i,k}).^2)* h);
    end
  
end

fig0 = figure;
fig1 = figure;
color = hsv(length(m));

for k = 1: length(m)
    figure(fig0)
    strdd = sprintf('Diamond Differenc2 S{%i}', m(k));
    strsc = sprintf('Step Characteristics S{%i}', m(k));
    
    plot(1 ./n, err0dd{k}, '.--','DisplayName', strdd, 'LineWidth', 2,...
        'Color', color(k, :), 'MarkerSize', 15) 
    hold on
    plot( 1 ./ n, err0sc{k}, '-','DisplayName', strsc, 'LineWidth', 2,...
        'Color', color(k, :),'Marker','.', 'MarkerSize', 15)
   
    
    figure(fig1)
    plot(1 ./n, err1dd{k}, '.--','DisplayName', strdd, 'LineWidth', 2,...
        'Color', color(k, :), 'MarkerSize', 15) 
    hold on
    plot( 1 ./ n, err1sc{k}, '.-', 'DisplayName', strsc, 'LineWidth', 2,...
        'Color', color(k, :), 'MarkerSize', 15)
    
end

figure(fig0)
title('L^2 error between Diffusion and Spatial Discretizations for \phi_0', 'FontSize', 15)
legend

figure(fig1)
title('L^2 error between Diffusion and Spatial Discretizations for \phi_1', 'FontSize', 15)
legend

end