function [phi0, phi1, psil] = step_characteristics(edges, slab, psi_l,...
    psi_r, ang_flag, k, tol)

%Given a slab discretized by 'edges' and decribed 'mat' for each sub
%-domain we perform step characteristics transport sweep using the boundary 
%conditions psi0_r, psi1_l and the angular discretization given by ang_flag
%and k until solutions are within a given tolerance tol
if length(edges) ~= (length(slab) + 1)
    error('Slab is improperly defined')
end

len = length(slab);
[Oz, w] = angles(ang_flag, k); %angular discretization
lenOz = length(Oz); 
psil = zeros(lenOz, len + 1); %for forward sweeps
psir = zeros(lenOz, len + 1); %for backward sweeps
psil(:, 1) = psi_l; %fill in left-side boundary
psir(:, len + 1) = psi_r; %fill in right-side boundary

%initial guesses for moment of scalar flux
phi0 = rand(1, len); 
phi1 = rand(1, len); 

%itermediary vectors
phi0_ang = zeros(lenOz, len);
phi1_ang = zeros(lenOz, len); 

err = sqrt(sum(phi0.^2) + sum(phi1.^2)); %error in calculated phis

max_iter = 3e03; %break while loop after max_iter iterations 
iter = 0; %iteration count
while err > tol && iter <= max_iter
    %old values to compare to new values
    
    for j = 1 : lenOz
        if Oz(j) > 0 %forward sweep 
            for k = 1: len
                obj = slab{k};
                
                %calculate outgoing flux and modified source
                [psil(j, k+1), Q] = obj.step_char(Oz(j), psil(j, k),...
                    phi0(k), phi1(k));
                
                %calculate PHI term for given value of Oz
                [phi0_ang(j, k), phi1_ang(j, k)] = obj.phi_maker(Oz(j),...
                    w(j), Q, psil(j,k), psil(j, k+1));
                
            end
                  
         elseif Oz(j) < 0 %backward sweep
             for k = 1: len
                 i = len + 1 - k; %iteration count for backward sweep
                 obj = slab{i};
                 
                 %calculate outgoing flux and modified source
                 [psir(j, i), Q] = obj.step_char(Oz(j), psir(j, i+1),...
                     phi0(i), phi1(i));
                 
                 %calculate PHI term for given value of Oz
                 [phi0_ang(j, i), phi1_ang(j, i)] = obj.phi_maker(Oz(j),...
                     w(j), Q, psir(j, i+1), psir(j, i));
                
             end
        end
    end
    
    oldphi0 = phi0;
    oldphi1 = phi1;
    phi0 = sum(phi0_ang); %sum over angles 
    phi1 = sum(phi1_ang); %to update PHIs
    
    %update error
    err = sqrt( sum((oldphi0 - phi0).^2) + sum((oldphi1 - phi1).^2) );
    
    iter = iter + 1; %update iteration count 
    
    if iter == max_iter
        disp(['Final error is' num2str(err)])
        error('Maximum number of iterations achieved')
    end
           
end



end 