function [] = gaussian_asy_main()
    % constants
    H12 = [1,0,0,0;
            0,-1,0,0;
    		0,0,-1,0;
    		0,0,0,1];
       
    H1 = [0,0,1,0;
    	0,0,0,1;
    	1,0,0,0;
    	0,1,0,0];
    
    H2 = [0,1,0,0;
    	1,0,0,0;
    	0,0,0,1;
    	0,0,1,0];

    gz = pi/4;
    dt = 1/500;

    t = datetime;
    tdy = sprintf('%i%i%i',t.Year, t.Month, t.Day);
    
    time_entropy = csvread('entropy_time.csv');

    % variables
    percentage = 0.5;
    terms = 3;
    rn = 50;
    target_ent = log(2) - percentage*log(2);
    
    ind = (time_entropy(:,1) == percentage);
    non_local_timing = time_entropy(ind,3)-dt;
    disp(non_local_timing);
    
    filename = sprintf('%s_landscape_gaussian_asy.csv', tdy);
    currentfolder = pwd;
    subfolder = fullfile(currentfolder,'Gaussian');
    filename = fullfile(subfolder,filename);
    tol = 5e-4;

    if exist(filename,'file') == 0
        masterdf = ones(1,terms + 2 + 1)*tol;
        csvwrite(filename,masterdf);
    else
        dlmwrite(filename,ones(1,terms + 2 + 1)*tol,'-append');
    end
    
    for i = 1:1
        fprintf('run no: %i \n', i);
        th_init = rand(2*terms+2+1,2)*2*pi;
        amplitude = randi([-rn,rn],2*terms+2+1,2);
        mean = rand(2*terms+2+1,2)*2;
        sd = rand(2*terms+2+1,2)*rn;

        param = [th_init, amplitude, mean, sd];

        [local_timing,best_param] = gaussian_asy_simplex(param,target_ent,gz,...
            rn,H1,H2,H12,dt,tol,filename);

%         [non_local_timing, ~] = entropy(best_param,target_ent,dt);

        if local_timing < non_local_timing
            dlmwrite(fullfile(subfolder,sprintf('%s_results_gaussian_asy.csv',tdy)),...
                [local_timing, non_local_timing, target_ent, tol, best_param, ...
                mod(best_param(1:2)/pi,2)], '-append','precision',10);
        end            
    end
    fprintf('start time: %s \n',char(t));
    t = datetime;
    fprintf('end time: %s \n',char(t));
end