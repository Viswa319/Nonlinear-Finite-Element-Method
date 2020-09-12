%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
%% In element routine funtion calculates elemental stiffness tensor, internal force and strain 
function [Ke,stress,strain,Fint_e,stress_ov] = element_routine(Q,relem,uelem,delta_uelem,prev_stress_ov)
           %% Inputs: 
           %%             relem :- list of size 2x1 with outer and inner element radius
           %%             uelem :- list of size 2x1 with element displacement 
           %%             delta_uelem :- list of size 2x1 with element delta displacement
           %%             prev_stress_ov :- list of size 2x1 with element over stress 
           %% Outputs:
           %%             Ke :- element stiffness tensor of size 2x2
           %%             stress :- stress of element called from material routine of size 2x1
           %%             strain :- strain of element of size 2x1
           [E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters();
            B = Belem(relem,0); % called from Belem function below
            strain = B* uelem; % getting element strain from element displacement
            delta_strain = B*delta_uelem; % getting element delta strain from element delta displacement
            [stress,C,stress_ov] = material_routine(Q,strain,delta_strain,prev_stress_ov); % calling material routine
            J = Jacobian(relem,Guass_point); % calling from Jacobian function below
             Fint_e = weights*B'*stress* (shape_function(Guass_point)'*relem)*J; % calculating internal element force of size 2x1
            Ke = weights*B'*C*B* (shape_function(Guass_point)'*relem)*J; % elemental stiffness tensor of size 2x2     
end
function N = shape_function(s)
           %% Shape function, lagrange polynomial of two elements in one dimensional
           %% Input: s :- elemental coordinate
           %% Output: returns a shape function   
            N = (1/2)*[1-s;1+s];         
  end
function dN = der_shape_function(s)
           %% derivative of shape function with respect to elemental coordinates
           %% Input: s :- elemental coordinate
           %% Output: returns a derivative of shape function w.r.t elemental coordinates
            dN = (1/2)*[-1;1];   
 end
function J = Jacobian(relem,s)
            %% Jacobian, derivative of global coordinates to elemental coordinates
            %% Input: relem :- list of size 2x1 with outer and inner element radius
            %% Output: Jacobian   
            dN = der_shape_function(s);
            J = dN'*relem;
 end
function B = Belem(relem,s)
            %% Parameter which relates strain and displacement
            %% Input: relem :- list of size 2x1 with outer and inner element radius and elemental coordinates
            %% Output: 2x2 matrix which relates strain and displacement like strain = B*displacement
             helem = relem(2)-relem(1);
             B =[-1/helem , 1/helem; (1-s)/(relem(1)*(1-s)+relem(2)*(1+s)) , (1+s)/(relem(1)*(1-s)+relem(2)*(1+s))]; 
end