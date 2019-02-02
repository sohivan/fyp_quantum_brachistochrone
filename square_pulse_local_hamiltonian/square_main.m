function [] = square_main()
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
    
    filename = sprintf('%s_landscape_square_s1xs2z.csv', tdy);
    currentfolder = pwd;
    subfolder = fullfile(currentfolder,'Square');
    filename = fullfile(subfolder,filename);
    tol = 5e-4;

    if exist(filename,'file') == 0
        masterdf = ones(1,terms + 4 + 4)*tol;
        csvwrite(filename,masterdf);
    else
        dlmwrite(filename,ones(1,terms + 4 + 4)*tol,'-append');
    end
    
    for i = 1:10
        fprintf('run no: %i \n', i);
        th_init = rand(terms+5,2)*2*pi;
        ps_init = zeros(terms+5,2);
        coef = randi([0,rn],terms+5,terms);

        param = [th_init, ps_init, coef];
        param(:,6) = -1*param(:,6);

        [local_timing,best_param] = square_simplex(param,target_ent,...
            gz,rn,H12,HXX,dt,terms,tol,filename);

%         [non_local_timing, ~] = entropy(best_param,target_ent,dt);

        if local_timing < non_local_timing
            dlmwrite(fullfile(subfolder,sprintf('%s_results_square_s1xs2z.csv',tdy)),...
                [local_timing, non_local_timing, target_ent, tol, best_param,...
                mod(best_param(1:4)/pi,2)], '-append','precision',10);
        end            
    end
    fprintf('start time: %s \n',char(t));
    t = datetime;
    fprintf('end time: %s \n',char(t));
end