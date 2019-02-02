function [t_end,max_ent] = square_asy_evolver(coef,ent_targ,gz,H1,H2,H12,dt)
% 	This function evolves the initial state and outputs a final state
% 	hardcode H12.HXX,gz,w for now
% 	t_holder as a vector, linspace()
%     gz*H12
%     add function entropy full into the code 
    
    max_ent = -Inf;

    k_init = productstate(coef(1),coef(2),0,0);
    
    iter = 2/dt;
    
    pos_1 = coef(3);
    neg_1 = coef(5);
    freq_1 = coef(7);
    
    pulse_1 = square(linspace(0,freq_1*pi,iter));
    pulse_1 = (pulse_1*pos_1).*(pulse_1>0) + (pulse_1*neg_1).*(pulse_1<=0);
    
    pos_2 = coef(4);
    neg_2 = coef(6);
    freq_2 = coef(8);
    
    pulse_2 = square(linspace(0,freq_2*pi,iter));
    pulse_2 = (pulse_2*pos_2).*(pulse_2>0) + (pulse_2*neg_2).*(pulse_2<=0);
    
    t_holder = 0;
       
    for j = 1:iter

		fn_1 = pulse_1(j);
        fn_2 = pulse_2(j);
		Htot = gz*H12 + fn_1*H1 +fn_2*H2;
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
            % disp('helper');
            % disp(ent_targ);
            % disp(ent_holder);
            % disp(t_holder);            
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