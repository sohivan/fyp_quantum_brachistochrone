function [t_end,max_ent] = square_evolver(coef,ent_targ,gz,H12,HXX,dt)
% 	This function evolves the initial state and outputs a final state
% 	hardcode H12.HXX,gz,w for now
% 	t_holder as a vector, linspace()
%     gz*H12
%     add function entropy full into the code 
    
    max_ent = -Inf;

    k_init = productstate(coef(1),coef(2),coef(3),coef(4));
    
    iter = 2/dt;
    
    pos = coef(5);
    neg = coef(6);
    freq = coef(7);
    
    pulse = square(linspace(0,freq*pi,iter));
    pulse = (pulse*pos).*(pulse>0) + (pulse*neg).*(pulse<=0);
    
    t_holder = 0;
    
    for j = 1:iter

		fn = pulse(j);
		Htot = gz*H12 + fn*HXX;
		Utot = expm(-1i*Htot*dt);
		k_init = Utot*k_init;

        density_matrix12 = (k_init)*(k_init');
        
        a = density_matrix12(1,1) + density_matrix12(2,2);
        b = density_matrix12(1,3) + density_matrix12(2,4);
        c = density_matrix12(3,1) + density_matrix12(4,2);
        d = density_matrix12(3,3) + density_matrix12(4,4);

        lambda_1 = 0.5*(a + d  - sqrt((a^2) + (4*b*c) - (2*a*d) + (d^2)));
        lambda_2 = 0.5*(a + d  + sqrt((a^2) + (4*b*c) - (2*a*d) + (d^2)));

        ent_holder = real(- (lambda_1*log(lambda_1) + lambda_2*log(lambda_2)));
        
        if ent_holder > max_ent
            max_ent = ent_holder;
        end
        
        if ent_targ < ent_holder
            break
        end

	    t_holder = t_holder + dt;
    end

    if isequal(j, iter)
        t_end = 1/max_ent;
    else
        t_end = t_holder;
    end
end