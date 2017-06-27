function driver_generate_assemblyFile()
    % Path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Path to the project directory
    directory_project = '../';
    
    
    addpath(directory_IGARoutines);
    addpath(directory_project);
    
    numRefinements = [0; 1; 2];
    
    for i = 1 : size(numRefinements, 1)
        generate_assemblyFile(numRefinements(i));
    end
    
    
    rmpath(directory_IGARoutines);
    rmpath(directory_project);
end
