function driver_3d_refine1()
    % Set path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Set path to the model directory
    directory_model = '../../../../models/j2_plasticity/';
    
    % Add paths
    addpath(directory_IGARoutines);
    addpath(directory_model);
    
    
    % Parameters for simulation
    numRefinements = 1;
    
    
    for i = 1 : size(numRefinements, 1)
        % Set path to the assembly directory
        path_to_assembly_directory = sprintf('../assembly_files/numRefinements%d/', numRefinements(i));
        
        % Set path to the results directory
        path_to_results_directory = sprintf('../results/numRefinements%d/', numRefinements(i));
        
        % Create the results directory if it does not exist
        if ~exist(path_to_results_directory, 'dir')
            mkdir(path_to_results_directory);
        end
        
        
        % Run the constitutive model
        model_3d(path_to_assembly_directory, path_to_results_directory, 0);
    end
    
    
    % Remove paths
    rmpath(directory_IGARoutines);
    rmpath(directory_model);
end
