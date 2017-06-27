function driver_postprocess_part4()
    % Set path to the IGA routines directory
    directory_IGARoutines = '../../../../routines_for_iga/';
    
    % Set path to the post-process routines directory
    directory_postprocessRoutines = '../';
    
    % Add paths
    addpath(directory_IGARoutines);
    addpath(directory_postprocessRoutines);
    
    
    % Run the post-processing routine
    postprocess_part4();
    
    
    % Remove paths
    rmpath(directory_IGARoutines);
    rmpath(directory_postprocessRoutines);
end
