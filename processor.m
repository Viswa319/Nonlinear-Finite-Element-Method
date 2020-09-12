%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
%Input parameters
function [u,ur_b,stress,t] = processor(Q,rnodes)
            [E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters();
            t =  []; % Initializing empty list for time steps
            t(1) = delta_t; % 
            kmax = 40;
            % Assignment matrix
            A = cell(1,nelem);
            for j = 1:nelem
                        A{j} = zeros(2,nelem+1);
            end
            for i = 1:nelem
                        A{i}(1,i) = 1;
                        A{i}(2,i+1) = 1;
            end
            function inf_norm= inf_norm(x)
                        %% function which gives infinity norm of a vector
                        %% 
                        inf_norm = max(abs(x));
            end
            stress = zeros(nelem+1,2); %initializing stress values as zeros of size nelem x 2 where each row has radial and hoop stress for respective element
            strain = zeros(nelem+1,2); %initializing strain values as zeros of size nelem x 2 where each row has radial and hoop strain for respective element
            stress_ov = zeros(nelem+1,2); %initializing stress values as zeros of size nelem x 2 where each row has radial and hoop over stress for respective element
            delta_u = zeros(nelem+1,1); %initializing delta displacement values as zeros of size nelem+1 x 1 
            u = zeros(nelem+1,1); %initializing displacement values as zeros of size nelem+1 x 1
            ur_b = [];
            ur_b(1) = 0;
            m = 2;
            for T1= delta_t:delta_t:tf
                        t(m)  = t(m-1)+delta_t;
                        if T1 <= tL
                                    p = Pmax*T1*(1/tL);
                        else
                                    p = Pmax;
                        end  
                        k = 1;
                        while  true
                                    F_int = zeros(nelem+1,1); 
                                    K = zeros(nelem+1,nelem+1); % creating a zero matrix for global stiffness matrix
                                    Ke = cell(1,nelem);
                                    Fext = zeros(nelem+1,1); % creating zero vector for external force
                                    Fext(1) = p*a; % Initializing first term in the Fext to p*a, where p increases with time
                                    for i = 1:nelem
                                                relem = rnodes(i:i+1); % nodal radius of corresponding element
                                                uelem = u(i:i+1); % dispalcement for corresponding element
                                                delta_uelem = delta_u(i:i+1); % delta displacement element
                                                % calling from element routine
                                                [Ke{i},stress_elem,strain_elem,Fint_elem,stress_ov_elem] = element_routine(Q,relem,uelem,delta_uelem,stress_ov(i,1:2)'); 
                                                K = K + (A{i}'*Ke{i}*A{i}); % Global stiffness tensor
                                                F_int=F_int + (A{i}'*Fint_elem); % Global internal force
                                                stress(i,1:2) = stress_elem(1:2); % storing stress values
                                                strain(i,1:2) = strain_elem(1:2); % storing strain values
                                                stress_ov(i,1:2) = stress_ov_elem(1:2); % updating over stress values
                                    end
                                    delta_u = linsolve(K,Fext-F_int); % delta displacement from Newton-Raphson method
                                    u = u+delta_u; % displacement
                                    ur_b(m) = u(end);
                                    if  inf_norm((Fext-F_int)) <0.005*(inf_norm(F_int))||inf_norm(delta_u)<0.005*(inf_norm(u))
                                                %fprintf("converged for iteration:");
                                                k;
                                                break
                                    end
                                    k = k+1;
                        end               
                        m = m+1;
            end 
            strain
end