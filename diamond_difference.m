function [phi0, phi1, psi] = diamond_difference(edges, mat, psi_l,...
    psi_r, ang_flag, k, tol)
%Given a slab discretized by 'edges' and decribed 'mat' for each sub
%-domain we perform diamond difference transport sweep using the boundary 
%conditions psi0_r, psi1_l and the angular discretization given by ang_flag
%and k until solutions are within a given tolerance tol
if length(edges) ~= (length(mat) + 1)
    error('Slab is improperly defined')
end

len = length(mat);
Oz = angles(ang_flag,k); %angular discretization
psi = zeros(length(Oz), len + 1); %output angular flux
%fill in boundary conditions of psi
psi(:, 1) = psi_l;
psi(:, len) = psi_r;

%initial guesses for moment of scalar flux
phi0 = rand(1, len); 
phi1 = rand(1, len); 

%itermediary vectors
phi0_ang = zeros(length(Oz), len);
phi1_ang = zeros(length(Oz), len); 

err = sqrt(sum(phi0.^2) + sum(phi1.^2)); %error in calculated phis

max_iter = 1e01; %break while loop after max_iter iterations 
iter = 0; %iteration count
while err > tol && iter <= max_iter 
   
    
    for j = 1 : length(Oz)
        if Oz(j) > 0 %forward sweep 
            for k = 1: len
                obj = mat{k};
                Q = obj.Q0;
                [phi0_ang(j, k), phi1_ang(j, k), Q, psi(j, k+1)] = ...
                    obj.diamond_difference(Oz(j), psi(j, k), phi0(k),...
                    phi1(k), Q);
                mat{k}.Q0 = Q; %update material source
            end
                  
         elseif Oz(j) < 0 %backward sweep
             for k = 1: len
                 i = len + 1 - k; %iteration count for backward sweep
                 obj = mat{i};
                 Q = obj.Q0;
                  [phi0_ang(j, i), phi1_ang(j, i), Q, psi(j, i)] = ...
                      obj.diamond_difference(Oz(j), psi(j, i+1),...
                      phi0(i), phi1(i), Q); 
                  mat{i}.Q0 = Q; %update material source
             end
        end
    end
    oldphi0 = phi0;
    oldphi1 = phi1;
    phi0 = sum(phi0_ang); %sum over angles (WE NEED THE WEIGHTS!!!)
    phi1 = sum(phi1_ang); %to update phis
    
    %update error
    err = sqrt( sum((oldphi0 - phi0).^2) + sum((oldphi1 - phi1).^2) );
    
    iter = iter + 1; %update iteration count 
    
    if iter == max_iter
        error('Maximum number of iterations achieved')
    end
           
end
     
end
 