%% Sai Viswanadha Sastry, Upadhyayula
%% 65130
%% Nonlinear Finite Element Method Assignment
%% 01/07/2020
%% lecturer in charge: Dr. Geralf Hütter
%% Generate list of position of nodes according to a geometric series
function rnodes = meshGenerator(a,b,nelem)
           % Function gives node points as a list
           %% Input: 
           %%        a :- inner radius
          %%         b :- outer radius
          %%         nelem :- number of elements
           % Output: returns a vector with nodal points for given number of elements
           if  nelem == 1
                      rnodes = [a;b];
            else
                        meshrefinementfactor=2; %ratio of element sizes at outer and inner radius
                        %ratio between element sizes of subsequent elements for a geometric series
                         q=meshrefinementfactor^(1./(nelem-1));
                        %size of first interval
                        dr=(b-a)*(1-q)/(1-meshrefinementfactor*q);
                        rnode=a;
                        rnodes=[a];
                        %loop over all elements
                        for i=1:nelem
                                    rnode=rnode+dr;
                                    rnodes=[rnodes;rnode];
                                    dr=dr*q;
                        end
             end
end