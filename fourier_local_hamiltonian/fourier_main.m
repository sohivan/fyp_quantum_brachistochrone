function [] = fourier_main()
    % constants
    H12 = [1,0,0,0;
            0,-1,0,0;
    		0,0,-1,0;
    		0,0,0,1];

%     H12 = [0,0,1,0;
%            0,0,0,-1;
%            1,0,0,0;
%            0,-1,0,0];

    HXX = [0,1,1,0;
    		1,0,0,1;
    		1,0,0,1;
    		0,1,1,0];

    gz = pi/4;
    w = pi/4;
    dt = 1/500;

    time_entropy = csvread('entropy_time.csv');

    t = datetime;
    tdy = sprintf('%i%i%i',t.Year, t.Month, t.Day);

    % variables
    percentage = 0.5;
    terms = 3;
    rn = 10;
    target_ent = log(2) - percentage*log(2);
    
    ind = (time_entropy(:,1) == percentage);
    non_local_timing = time_entropy(ind,3)-dt;
    
    filename = sprintf('%s_landscape_fourier_%i_s1xs2z.csv', tdy,terms);
    currentfolder = pwd;
    subfolder = fullfile(currentfolder,'Fourier');
    filename = fullfile(subfolder,filename);
    tol = 5e-4;

    if exist(filename,'file') == 0
        masterdf = ones(1,2*terms + 4 + 2)*tol;
        csvwrite(filename,masterdf);
    else
        dlmwrite(filename,ones(1,2*terms + 4 + 2)*tol,'-append');
    end
    
    for i = 1:500
        fprintf('run no: %i \n', i);
        th_init = rand(2*terms+5,2)*2*pi;
        ps_init = zeros(2*terms+5,2);
        coef = randi([-rn,rn],2*terms+5,2*terms);

        param = [th_init, ps_init, coef];
        param(:,6) = 0;

        [local_timing,best_param] = fourier_simplex(param,target_ent,w,gz,rn,H12,HXX,dt,tol,terms,filename);

%         [non_local_timing, ~] = entropy(best_param,target_ent,dt);

        if local_timing < non_local_timing
            dlmwrite(fullfile(subfolder,sprintf('%s_results_fourier_%i_s1xs2z.csv',tdy,terms))...
                , [local_timing, non_local_timing, target_ent, tol, best_param, ...
                mod(best_param(1:4)/pi,2)], '-append','precision',10);
        end            
    end
    fprintf('start time: %s \n',char(t));
    t = datetime;
    fprintf('end time: %s \n',char(t));
end