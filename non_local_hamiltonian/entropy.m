function [t_end, ent_max] = entropy(coefficient,ent_target,dt)
% this function calculates the entanglement entropy of a specific state
% it returns the entanglement entropy that is above the ent_target and the
% time needed for the particular state to reach that entanglement
%
% For instances when the particular state doesn't reach the target
% entanglement, what is the time I should return

    th1 = coefficient(1);
    th2 = coefficient(2);

    gz = pi/4;
    t_holder = 0;
    
    ent_max = -Inf;

    maxiter = 2/dt;

    for i = 1:maxiter

        tau = t_holder*gz;

        sqr = (cos(th1)^2) + (sin(th1)^2)*((cos(2*tau)^2) + (cos(th2)^2)*(sin(2*tau)^2));
        sqr = sqrt(sqr);

        ent_cur = (-0.5-0.5*(sqr))*log(0.5*(1+sqr))-(0.5-0.5*(sqr))*log(0.5*(1-sqr));

        if ent_max < ent_cur
            ent_max = ent_cur;
        end
        
        if ent_cur > ent_target
            break
        end
        
        t_holder = dt + t_holder;
        
    end
    
    if isequal(i,maxiter)
        t_end = 1/ent_max;
    else
        t_end = t_holder;
    end
    
end