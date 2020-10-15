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
        phi0 
        phi1
        phi2
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
               obj.nu = 4;
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
        
        
        function [psi1, phi0new, phi1new, phi2new, Q_new] = ...
                step_characteristic(obj,Oz, psi0, phi0, phi1, phi2, Q)
            
            %Takes the direction of neutron flow, incoming angular flux and
            %0th to 2nd scalar fluxes and performs a step characteristics
            %transport sweep to determine the angular flux within material 
            % and update the 0th to 2nd angular moments of the angular flux
            
            Delta = abs(obj.right - obj.left);
            tau = obj.sig_t * Delta /abs(Oz); 
            
            psi1 = psi0 * exp(-tau ) + Q / obj.sig_t * (1 - exp( -tau )); 
            
            
            
            
        end
        
        
        function [psi1, phi0new, phi1new, phi2new, Q_new] = ...
                diamond_difference(obj,Oz, psi0, phi0, phi1, phi2)
            
            %Takes the direction of neutron flow and the incoming angular
            %flux and performs a diamond difference transport sweep to
            %determine angular flux within material and the 0th to 2nd 
            %angular moments of the angular flux
            
            Delta = (obj.right - obj.left);
            tau = obj.sig_t * Delta / abs(Oz);
         
        end
            
    end
    
    methods(Static) 
        function P2_Oz = P2(Oz)
            %second order lagrange polynomial
            P2_Oz = (3 * Oz.^2 - 1 ) / 2;
  
        end
        
        
    end
    
    
end


