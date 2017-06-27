function driver_generate_assemblyFile()
    % Path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Path to the project directory
    directory_project = '../';
    
    
    addpath(directory_IGARoutines);
    addpath(directory_project);
    
    
    orders = [2; 4; 6; 8];
    numRefinements = [0; 1; 2; 3; 4];
    
    for i = 1 : size(orders, 1)
        for j = 1 : orders(i) / 2
            for k = 1 : size(numRefinements, 1)
                generate_assemblyFile(orders(i), j, numRefinements(k));
            end
        end
    end
    
    
    rmpath(directory_IGARoutines);
    rmpath(directory_project);
end
