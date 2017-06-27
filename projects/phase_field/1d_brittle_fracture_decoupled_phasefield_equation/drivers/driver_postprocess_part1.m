function driver_postprocess_part1()
    % Set path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Set path to the post-process routines directory
    directory_postprocessRoutines = '../';
    
    % Add paths
    addpath(directory_IGARoutines);
    addpath(directory_postprocessRoutines);
    
    
    % Parameters for simulation
    orders = [2; 4; 6; 8];
    numRefinements = [0; 1; 2; 3; 4];
    
    
    for i = 1 : size(numRefinements, 1)
        for j = 1 : size(orders, 1)
            % Set path to the assembly directory
            path_to_assembly_directory = sprintf('../assembly_files/order%d/numRefinements%d/', orders(j), numRefinements(i));
            
            % Set path to the results directory
            path_to_results_directory = sprintf('../results/order%d/numRefinements%d/', orders(j), numRefinements(i));
            
            % Set path to the outputs directory
            path_to_outputs_directory = strcat(path_to_results_directory, 'plots/');
            
            % Create the outputs directory if it does not exist
            if ~exist(path_to_outputs_directory, 'dir')
                mkdir(path_to_outputs_directory);
            end
            
            
            % Run the post-processing routine
            postprocess_part1(path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory, numRefinements(i));
        end
    end
    
    
    % Remove paths
    rmpath(directory_IGARoutines);
    rmpath(directory_postprocessRoutines);
end
