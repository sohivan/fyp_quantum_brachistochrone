function [best_time,best_param] = square_asy_simplex(param,ent,gz,rn,H12,H1,H2,dt,terms,tol,filename)

	% parameters for simplex behaviour from internet
	ap = 1;
	bt = 2;
	gm = 0.5;
	dl = 0.5;

	vertices = size(param,1);
	timeholder = zeros(vertices,1);
	s = zeros(vertices,1);

	for i = 1:vertices
		coef_i = param(i,:);
		[timeholder(i),s(i)] = square_asy_evolver(coef_i,ent,gz,H1,H2,H12,dt);
	end

	[timeholder,ind] = sort(timeholder);
	param = param(ind,:);
	s = s(ind);

	for k = 1:1000
		
		% getting the different values for comparison
		time_worst = timeholder(end);
		time_best = timeholder(1);
		
		% to calculate centroid
		cnt = (sum(param,1) - param(end,:))/(vertices-1);
		
		% reflect
		cref = cnt + ap*(cnt - param(end,:));
		[time_refl,s_refl] = square_asy_evolver(cref,ent,gz,H1,H2,H12,dt);
		if time_best <= time_refl && time_refl < timeholder(end-1) %intermediate
			param(end,:) = cref;
			timeholder(end) = time_refl;
			s(end) = s_refl;
			
		% reflect best
		elseif time_refl < time_best
			   % expand
			   cexp = cnt + bt*(cref-cnt);
			   [time_exp,s_exp] = square_asy_evolver(cexp,ent,gz,H1,H2,H12,dt);
			   if time_exp < time_refl
				   % expand is better
				   param(end,:) = cexp;
				   timeholder(end) = time_exp;
				   s(end) = s_exp;
			   else
				   % reflect is better
				   param(end,:) = cref;
				   timeholder(end) = time_refl;
				   s(end) = s_refl;
			   end
			   
		% bad intermediate
		elseif timeholder(end-1) <= time_refl && time_refl < time_worst
				% contract outwards
				ccon = cnt + gm*(cref-cnt);
				[time_con,s_con] = square_asy_evolver(ccon,ent,gz,H1,H2,H12,dt);
				% contract better
				if time_con <= time_refl
				   param(end,:) = ccon;
				   s(end) = s_con;
				   timeholder(end) = time_con;
				% shrink
				else
					param(2:vertices,:) = param(1,:) + dl*(param(2:vertices,:)-param(1,:));
					for i = 1:vertices
						coef_i = param(i,:);
						[timeholder(i), s(i)] = square_asy_evolver(coef_i,ent,gz,H1,H2,H12,dt);
					end
				end
				
		% reflect worst
		elseif time_worst <= time_refl
				% contract inwards
				ccon = cnt - gm*(cref-cnt);
				[time_con,s_con] = square_asy_evolver(ccon,ent,gz,H1,H2,H12,dt);
				% contract better
				if time_con < time_worst
				   param(end,:) = ccon;
				   timeholder(end) = time_con;
				   s(end) = s_con;
				% shrink
				else
					param(2:vertices,:) = param(1,:) + dl*(param(2:vertices,:)-param(1,:));
					for i = 1:vertices
						coef_i = param(i,:);
						[timeholder(i), s(i)] = square_asy_evolver(coef_i,...
                            ent,gz,H1,H2,H12,dt);
					end                  
				end       
		end

		[timeholder,ind] = sort(timeholder);
		param = param(ind,:);
		s = s(ind);

		% finding the standard deviation of the vertices
		strd = mean(std(param));

		if strd < tol
			th_rand = rand(1,2)*2*pi;
			coef_rand = randi([1,rn],1,terms);
            coef_rand(5:6) = -1*coef_rand(5:6);
			param_rand = [th_rand, coef_rand];
			[time_rand, s_rand] = square_asy_evolver(param_rand,ent,gz,...
                H1,H2,H12,dt);

			if time_rand > timeholder(1)
			% if the random point doesn't result in a better timing, we terminate the search
				break
			else
			% if the random point does result in a better timing, we will use that new point as the best new point
				timeholder(1) = time_rand;
				param(1,:) = param_rand;
				s(1) = s_rand;
			end
		end

	end
	fprintf('no of iteration: %i \n', k);
	best_param = param(1,:);
	best_s = s(1);
	best_time = timeholder(1);

	store = [best_s, best_time, ent, tol, best_param, mod(best_param(1:2)/pi,2)];
	dlmwrite(filename,store,'-append','precision',10);
end
