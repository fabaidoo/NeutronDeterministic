function  SlabSystem
%Symmetric composite slab: Absorber (1cm), Reflector (10 cm), Isotropic (1
%cm), Scatterer (4 cm), Fuel (2 cm), Scatterer (4cm), Isotropic (1cm),
%Reflector (10 cm), Absorber (1 cm). Gauss-Legendre angular discretization 
%for n = 2, 6, 12 and 20. System spatially discretized with uniform number 
%of meshes.

%LOCATION OF EDGES OF MATERIAL in SYSTEM
mat_edges = [0 1 11 12 16 18 22 23 33 34];                                  %[0 1 1.5 4 4.5 5.5];
%MATERIALS IN SYSTEM
mat = {'absorber', 'reflector', 'isotropic', 'scatterer', 'fuel',...        %{'isotropic', 'absorber', 'fuel', 'absorber', 'isotropic'};
    'scatterer', 'isotropic', 'reflector', 'absorber'};
%EDGES OF MESH
h = .05; %mesh size. Chosen so meshes coincide with material edges
edges = mat_edges(1) : h : mat_edges(length(mat_edges));
z = edges(1: length(edges) - 1) + h / 2 ; %domain over which we plot

%Create SLAB
slab = cell(1, length(edges) - 1);
for j = 1: length(mat)
    %indices of material edges in slab
    i_l = find(edges == mat_edges(j)); 
    i_r = find(edges == mat_edges(j+1));
   for k = i_l : i_r - 1  
       slab{k} = material(mat{j}, edges(k), edges(k + 1));
   end
end

%ANGULAR DISCRETIZATION
flag = 'glq';
n = [2, 6, 12, 20];
%TOLERANCE
tol = 1e-4;

fig0 = figure;
fig1 = figure;
color = hsv(length(n)); %color map for plots

for i = 1: length(n)
    
    %VACUUM BOUNDARY CONDITIONS ON LEFT AND RIGHT
    psil = zeros(n(i) ,1);
    psir = zeros(n(i), 1); 
    
    [phi0sc, phi1sc] = step_characteristics(edges, slab, psil, psir,...
        flag, n(i), tol);
    
    [phi0dd, phi1dd] = diamond_difference(edges, slab, psil, psir,...
        flag, n(i), tol);
   
    strdd = sprintf('DD GLQ n = %i', n(i));
    strsc = sprintf('SC GLQ n = %i', n(i));
    figure(fig0)
    plot(z, phi0dd,'-', 'LineWidth', 2, ...                                  'Marker', '.', 'MarkerSize', 15,
        'Color', color(i, :), 'DisplayName', strdd)
    hold on
    plot(z, phi0sc,'--', 'LineWidth', 2,  ...                                'Marker','.','MarkerSize', 15,
        'Color', color(i, :), 'DisplayName', strsc)
   
    figure(fig1)
    plot(z, phi1dd,'-', 'LineWidth', 2, ...                                 'Marker', '.','MarkerSize', 15,
        'Color', color(i, :), 'DisplayName', strdd)
    hold on
    plot(z, phi1sc,'--', 'LineWidth', 2, ...                                'Marker','.', 'MarkerSize', 15,
        'Color', color(i, :), 'DisplayName', strsc)
end

str0 = sprintf(' for Symmetric slab (mesh = %.3f cm)', h);
str1 = sprintf(' for Symmetric slab (mesh = %.3f cm)', h);

figure(fig0)
plot(mat_edges, zeros(1, length(mat_edges)),  'k|-','MarkerSize', 35,...
    'DisplayName', 'Material Edges')
title(['\phi_0' str0],'FontSize', 18)
lgd0 = legend;
lgd0.FontSize = 14;
t=gca; t.XAxis.TickLength = [0 0];
%set(gca,'XTick',[])
xlabel('z (cm)')
xlim([mat_edges(1) mat_edges(length(mat_edges))])

figure(fig1)
plot(mat_edges, zeros(1, length(mat_edges)),  'k|-','MarkerSize', 30,...
    'DisplayName', 'Material Edges')
title(['\phi_1' str1], 'FontSize', 18 )
lgd1 = legend;
lgd1.FontSize = 14;
t=gca; t.XAxis.TickLength = [0 0];
xlabel('z (cm)')
xlim([mat_edges(1) mat_edges(length(mat_edges))])














    
    





end