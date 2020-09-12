%Input parameters
function [E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters()
            E = 200000; % Young's modulus in MPa
            nu = 0.2; % Poisson's ratio
            T = 1; % Time scale in sec
            a = 50; % Inner radius in mm
            b = 100; % Outer radius in mm
            Pmax = 140; % Maximum internal pressure to apply in MPa
            tL = 2; % Time from where pressure is kept constant in sec
            tf = 10; % Time until pressure is hold in sec
            nelem = 15; % Total number of elements
            delta_t = 0.1; % time step
            rnodes = meshGenerator(a,b,nelem);
            weights = 2;
            Guass_point = 0;
 end