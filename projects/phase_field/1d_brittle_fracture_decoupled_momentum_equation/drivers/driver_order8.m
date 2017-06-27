function driver_order8()
    % Set path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Set path to the model directory
    directory_model = '../../../../models/phase_field/1d_brittle_fracture_decoupled_momentum_equation/';
    
    % Add paths
    addpath(directory_IGARoutines);
    addpath(directory_model);
    
    
    % Parameters for simulation
    order = 8;
    numRefinements = [0; 1; 2; 3; 4];
    
    
    for i = 1 : order / 2
    for j = 1 : size(numRefinements, 1)
        % Set path to the assembly directory
        path_to_assembly_directory = sprintf('../assembly_files/order8_p%d/numRefinements%d/', i, numRefinements(j));
        
        % Set path to the results directory
        path_to_results_directory = sprintf('../results/order8_p%d/numRefinements%d/', i, numRefinements(j));
        
        % Create the results directory if it does not exist
        if ~exist(path_to_results_directory, 'dir')
            mkdir(path_to_results_directory);
        end
        
        
        % Run the constitutive model
        model_1d(order, path_to_assembly_directory, path_to_results_directory);
    end
    end
    
    
    % Remove paths
    rmpath(directory_IGARoutines);
    rmpath(directory_model);
end
