%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
function  [stress,C,stress_ov] = material_routine(Q,strain,delta_strain,prev_stress_ov)
           %% Function which returns stress and over stress. Over stress is calculated using Euler backward method
           %% Inputs: strain, delta strain, previous over stress
            %% Outputs: stress, tangential stiffness tensor, overstress 
           [E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters();
            Ce= (E/((1+nu)*(1-2*nu)))*[1-nu , nu;nu,1-nu];   % Elastic stiffness tensor
            Ct = Q*(1/(1+(delta_t/(2*T))))* [2/3, -1/3;-1/3,2/3] ;% Visco-elastic stiffness tensor
            C = Ce + Ct; % Tangential stiffness tensor 
            hydrostatic_strain = ((sum(delta_strain))/3)*[1;1]; % Hydrostatic strain  
            dev_strain = delta_strain - hydrostatic_strain; % Deviatoric strain
            stress_ov = 1/(1+(delta_t/(T)))*(prev_stress_ov+(Q*dev_strain)); % over stress using Euler backward method
            stress = Ce*strain + stress_ov; % Overall stress
 end