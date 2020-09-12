%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
clc;
clear all;
Problem = input("Linear :- 0, Non-linear :- 1 :");            
if Problem == 0
            Q = 0;
else
             Q = 100000; % Modulus which governs over stress evolution in MPa
end
[E,nu,T,a,b,Pmax,tL,tf,nelem,delta_t,rnodes,weights,Guass_point] =  Input_parameters();
[u,ur_b,stress,t] = processor(Q,rnodes); % Desired values called from processor function
u_anal = analytical(rnodes); % analytical displacement called from analytical function
% Plots 
figure(1)
plot(rnodes,u_anal, 'r -',rnodes,u,'b x');
axis([rnodes(1)-1  rnodes(end)+1]);
title('Radius vs Radial displacement');
xlabel('radius(mm)');
ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
legend('analytical solution ','numerical solution','location','NorthEast');
title('Radius vs Radial displacement');
figure(2);
subplot(2,1,1);
plot(rnodes(1:end),stress(:,1),' v -o');
title('Radius vs Radial stress','Interpreter', 'tex');
xlabel('radius r (mm)');
ylabel('{\sigma}_{r r}(MPa)','Interpreter', 'tex');
axis([rnodes(1)-1  rnodes(end)+1 ]);
subplot(2,1,2);
plot(rnodes(1:end-1),stress(1:end-1,2),'r-*','linewidth',1);
title('Radius vs Hoop stress','Interpreter', 'tex');
xlabel('radius r (mm)');
ylabel('{\sigma}_{\phi \phi} (MPa)','Interpreter', 'tex');
axis([rnodes(1)-1  rnodes(end)+1]);
figure(3)
plot(t,ur_b);
title('Time history of the widening of the pipe u_{r}(r=b,t)','Interpreter', 'tex');
xlabel('Time t (s)');
ylabel('Radial displacement u_{r} (mm)','Interpreter', 'tex');
axis([t(1)-1  t(end)+1]);
Table = [rnodes,u,u_anal,stress(:,1),stress(:,2)];