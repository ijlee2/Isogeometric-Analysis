function driver_postprocess_refine1()
    % Set path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Set path to the post-process routines directory
    directory_postprocessRoutines = '../';
    
    % Add paths
    addpath(directory_IGARoutines);
    addpath(directory_postprocessRoutines);
    
    
    % Parameters for simulation
    numRefinements = 1;
    
    
    for i = 1 : size(numRefinements, 1)
        % Set path to the assembly directory
        path_to_assembly_directory = sprintf('../assembly_files/numRefinements%d/', numRefinements(i));
        
        % Set path to the results directory
        path_to_results_directory = sprintf('../results/numRefinements%d/', numRefinements(i));
        
        % Set path to the outputs directory
        path_to_outputs_directory = sprintf('../results/numRefinements%d/plots/', numRefinements(i));
        
        % Create the outputs directory if it does not exist
        if ~exist(path_to_outputs_directory, 'dir')
            mkdir(path_to_outputs_directory);
        end
        
        
        for time_index = 1 : 1 : 1000
            % Run the post-processing routine
            postprocess(path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory, time_index);
        end
    end
    
    
    % Remove paths
    rmpath(directory_IGARoutines);
    rmpath(directory_postprocessRoutines);
end
