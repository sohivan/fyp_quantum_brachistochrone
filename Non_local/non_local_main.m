function [] = non_local_main()
    % constants
%     H12 = [1,0,0,0;
%            0,-1,0,0;
%            0,0,-1,0;
%            0,0,0,1];
% 
%     gz = pi/4;
    dt = 1/500;
    
    t = datetime;
    tdy = sprintf('%i%i%i',t.Year, t.Month, t.Day);

    terms = 0;

    filename = sprintf('%s_entropy_time.csv', tdy);
    currentfolder = pwd;
    filename = fullfile(currentfolder,'Non_local',filename);
    tol = 5e-4;
    
    if exist(filename,'file') == 0
        masterdf = ones(1,terms + 4 + 4)*tol;
        csvwrite(filename,masterdf);
    else
        dlmwrite(filename,ones(1,terms + 4 + 4)*tol,'-append');
    end
    
    percentage = linspace(0.01,0.99,99);
    maxiter = 1000;
    
    for j = percentage
        target_ent = log(2) - j*log(2);
        fprintf('run no: %i \n', j);
        
        storer = zeros(maxiter,1);
    
        for i = 1:maxiter
%             fprintf('run no: %i \n', i);
            th_init = rand(terms+5,2)*2*pi;
            ps_init = zeros(terms+5,2);

            param = [th_init, ps_init];

            [storer(i),~] = simplex_non_local(param,target_ent,dt,tol,...
                filename);     
        end
        
        min_time = min(storer);
        
    	store = [j, target_ent, min_time];
    	dlmwrite(filename,store,'-append','precision',10);

    end
            
    fprintf('start time: %s \n',char(t));
    t = datetime;
    fprintf('end time: %s \n',char(t));
end