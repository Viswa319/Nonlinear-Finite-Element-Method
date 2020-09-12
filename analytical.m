%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
% Analytical solution for elastic problem
function u = analytical(rnodes)
           [E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters();
           u = zeros(length(rnodes),1);
           for i = 1:length(rnodes)
                        u(i) = (1+nu)*(Pmax/E)*(a^2/(b^2-a^2))*((1-2*nu)*rnodes(i)+(b^2)/rnodes(i));
end