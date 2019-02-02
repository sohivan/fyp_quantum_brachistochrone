clear

particle_no = 50;
parameter_no = 2;
inertia_coef = 0.8; % usually between 0.8 to 1.2
cognitive_coef = 2; % usually close to 2
social_coef = 2; % usually close to 2
maxiter = 2;

% initialise initial position of particles
x_init = rand(particle_no,parameter_no);

% initialise the best position of each particle
% the arrangement of the particles cannot be changed
personal_best = x_init;

% store the quality of each position
f_out = zeros(particle_no,1);

% storing snapshot at each iteration for view
snapshot = zeros(maxiter,particle_no,parameter_no);

for i = 1:particle_no
	x_i = x_init(i,:);
	f_out(i) = camel(x_i(1),x_i(2));
end

% sort out the quality of each position to obtain the best position
[f_sorted,ind] = sort(f_out);

% best position among the swarm and the best quality of the landscape
swarm_best = x_init(ind(1),:);
landscape_best = f_out(ind(1));

% initialise the velocity of each particle
v_init = rand(particle_no,parameter_no);

for j = 1:maxiter
	% randomize the weightage of social and cognitive coefficient
	r1 = rand(particle_no,parameter_no);
	r2 = rand(particle_no,parameter_no);

	% this is to store the position of particles as a snapshot of time
	snapshot(j,:,:) = x_init;

	% find the new velocity of each particle
	disp(size(r1.*(personal_best - x_init)));
	disp(size(r2.*(swarm_best - x_init)));
	v_init = inertia_coef*v_init + cognitive_coef*r1.*(personal_best - x_init) + social_coef*r2.*(swarm_best - x_init);
    disp(size(v_init));
	% move the particles based on this new velocity to the new position, unit of time is 1
	x_init = x_init + v_init;
    disp(size(x_init));

	% iterate through all the particles
	for i = 1:particle_no
		x_i = x_init(i,:);
		% check if the new position is better than the previous position for each particle
		% if yes, then update; else don't update
		f_intm = camel(x_i(1),x_i(2));
		if f_intm < f_out(i)
			personal_best(i,:) = x_i;
			% if the new position is better than the previous one, there is a possibility that it is also better than the leader
			% this checks for that
			if f_intm < landscape_best
				swarm_best = x_i;
				landscape_best = f_intm;
			end
		end
		f_out(i) = f_intm;
	end

	% need to add a suitable stopping criteria
end

% write the points into a csv file
% snapshot1 = reshape(snapshot,maxiter,[]);
% csvwrite('try.csv',snapshot1);