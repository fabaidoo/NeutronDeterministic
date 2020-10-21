classdef material
    %Material class contains properties of the material neutrons move
    %through as well as the spatially discretized methods to compute the
    %angular/scalar flux. Stacks of material constitute a slab. To
    %spatially discretize a slab, use row vector of material class whose
    %boundaries coincide with discretization
 
    
    properties
        sig_t = 0; %total cross-section
        sig_s0 = 0; %zeroth moment of scattering cross-section
        sig_s1 = 0; %first moment of scattering cross-section
        sig_m = 0; % multiplication cross-section
        nu = 0; %number of particles produced in multiplication reaction
        
        mattype %identifies type of material. Useful for debugging
        
        Q0 = 0; %zeroth moment of source
        Q1 = 0; %first moment of source
        
        left_bnd = 0; %location of left boundary
        right_bnd = 1; %location of right boundary  
        
        %angular moments of scalar flux 
        %phi0 
        %phi1
    end
    
    methods
        function obj = material(mat_type,left, right)
            %Specify type of material and location of left and right
            %boundary. Outputs material class with specified properties.
            if nargin == 3
               obj.left_bnd = left;
               obj.right_bnd = right; 
            end
                   
            if strcmpi(mat_type, 'fuel') == 1
               obj.mattype = 'fuel';
               obj.sig_t = 1;
               obj.sig_m = .2 ;
               obj.nu = 2.5;
               obj.sig_s0 = .6;
               obj.sig_s1 = .1;
            elseif strcmpi(mat_type, 'reflector') == 1
                obj.mattype = 'reflector';
                obj.sig_t = 2;
                obj.sig_s0 = 1.8;
                obj.sig_s1 = -1;
            elseif strcmpi(mat_type, 'scatterer') == 1 
                obj.mattype = 'scatterer';
                obj.sig_t = 2;
                obj.sig_s0 = 1.9;
            elseif strcmpi(mat_type,'absorber') == 1
                obj.mattype = 'absorber';
                obj.sig_t = 10;
                obj.sig_s0 = 2;
                obj.sig_s1 = 2;
            elseif strcmpi(mat_type, 'isotropic') == 1
                obj.mattype = 'isotropic';
                obj.sig_t = .100;
                obj.Q0 = 1 ;
            elseif strcmpi(mat_type, 'anisotropic') == 1
                obj.mattype = 'anisotropic';
                obj.sig_t = .1;
                obj.Q0 = 1;
                obj.Q1 = 1;
            elseif strcmpi(mat_type, 'pure absorber') == 1
                obj.mattype = 'pure absorber';
                obj.sig_t = 10;
            else
                error("Material type not available")
            end     
        end
        
        function [psi_out, Qnew] = diamond_diff(obj, Oz, psi_in, phi0, phi1)
            %takes incoming angular flux and previous 0th and 1st angular 
            %moments of average and calculates outgoing angular flux via
            %diamond difference method
            
            Delta = abs(obj.right_bnd - obj.left_bnd); %material width
            tau = obj.sig_t * Delta ./ abs(Oz); 
            
            %total source term
            Qnew = obj.Q0 + 0.5 * obj.sig_s0 * phi0 + 1.5 * obj.sig_s1...
                * phi1 * Oz + 0.5 * obj.nu * obj.sig_m * phi0;% / (4 * pi);
            
             %compute the exiting angular flux
            psi_out = psi_in .* (2 - tau) ./ (2 + tau) + ... 
                (Qnew ./ obj.sig_t) .* (1 - (2 - tau) ./ (2 + tau)); 
        end
    
        function [psi_out, Qnew] = step_char(obj, Oz, psi_in, phi0, phi1)
            %takes incoming angular flux and previous 0th and 1st angular 
            %moments of average and calculates outgoing angular flux via
            %step characteristics method
            
            Delta = abs(obj.right_bnd - obj.left_bnd); %material width
            tau = obj.sig_t * Delta ./ abs(Oz); 
            
            %total source term
            Qnew = obj.Q0 + 0.5 * obj.sig_s0 * phi0 + 1.5 * obj.sig_s1...
                * phi1 * Oz + 0.5 * obj.nu * obj.sig_m * phi0;% / (4 * pi);
            
            %compute the exiting angular flux
            psi_out = psi_in .* exp(-tau) + (Qnew ./ obj.sig_t) .* ...
                (1 - exp(-tau)); 
        end
        
        function [phi0_out, phi1_out] = phi_maker(obj, Oz, w, Q, psi0,...
                psi1)
            %Provides terms for new angular moments of phi for material. 
            %Takes angle Oz and its weight, incoming and outgoing psis
            %and modified source
            
            Delta = abs(obj.right_bnd - obj.left_bnd); %material width
            tau = obj.sig_t * Delta ./ abs(Oz); %has same shape as Oz 
            
            phi0_out =  (Q / obj.sig_t + (psi0 - psi1) / tau ) * w;
            
            phi1_out = (Q / obj.sig_t + (psi0 - psi1) / tau ) * Oz...
                * w;
        end
            
    end
      
    
end


