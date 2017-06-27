function generate_bcFile(order, numRefinements)
    % Set the path to the assembly files directory
    path_to_assembly_directory = sprintf('../assembly_files/order%d/numRefinements%d/', order, numRefinements);
    
    % File paths
    file_assembly_global = strcat(path_to_assembly_directory, 'file_assembly_global');
    file_assembly_patch  = strcat(path_to_assembly_directory, 'file_assembly_patch');
    
    % Load the global assembly file
    load(file_assembly_global, ...
         'numDOFs'       , ...
         'numDOFsPerNode', ...
         'GN_array');
    
    % Load patch 1
    load(sprintf('%s%d', file_assembly_patch, 1), ...
         'numNodesPerElement1', ...
         'numElements1'       , ...
         'IEN_array');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Set BCs for the phase field
    %   
    %   We let the phase field equal to 1 at the two ends of the bar by
    %   setting the first and last coefficients to 1. The derivative(s)
    %   of the phase field can be made equal to 0 by setting the next
    %   coefficients equal to 1 as well.
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Initialize the BC arrays
    BCs_displacement = [];
    BCs_force = [];
    
    
    % Degrees of freedom and BC values at x = -bar_length
    dofIndex = 1;
    bcValue = 1;
    
    elementIndex = 1;
    
    nodeIndex = find_node_index((1 : (order / 2))', ...
                                [], ...
                                [], ...
                                numNodesPerElement1, ...
                                [], ...
                                [], ...
                                IEN_array(:, elementIndex));
    
    temp = ones(order / 2, 1);
    
    BCs_displacement = [BCs_displacement; ...
                        [GN_array(nodeIndex), dofIndex * temp, bcValue * temp]];
    
    
    % Degrees of freedom and BC values at x = +bar_length
    dofIndex = 1;
    bcValue = 1;
    
    elementIndex = numElements1;
    
    nodeIndex = find_node_index(((numNodesPerElement1 - order / 2 + 1) : numNodesPerElement1)', ...
                                [], ...
                                [], ...
                                numNodesPerElement1, ...
                                [], ...
                                [], ...
                                IEN_array(:, elementIndex));
    
    BCs_displacement = [BCs_displacement; ...
                        [GN_array(nodeIndex), dofIndex * temp, bcValue * temp]];
    
    
    %----------------------------------------------------------------------
    %  Build the ID, LM, and BC arrays
    %----------------------------------------------------------------------
    ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array);
    LM_array2 = build_lm_array(IEN_array, ID_array, GN_array, 0);
    
    [BCU_array2, BCF_array2] = build_bc_array(BCs_displacement, BCs_force, ID_array);
    
    % Number of unknown degrees of freedom
    numUnknownDOFs2 = numDOFs - size(BCU_array2, 1);
    
    % Indices for the unknown coefficients
    index_u2 = (1 : numUnknownDOFs2)';
    index_f2 = ((numUnknownDOFs2 + 1) : numDOFs)';
    
    
    clear BCs_displacement BCs_force ID_array;
    
    
    
    %----------------------------------------------------------------------
    %  Save the BCs
    %----------------------------------------------------------------------
    save(sprintf('%sfile_bc', path_to_assembly_directory), ...
                 'LM_array2', ...
                 'BCU_array2', ...
                 'BCF_array2', ...
                 'numUnknownDOFs2', ...
                 'index_u2', ...
                 'index_f2', ...
                 '-v7.3');
end
