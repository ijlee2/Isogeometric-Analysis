function generate_bcFile(order, degree, numRefinements)
    % Set the path to the assembly files directory
    path_to_assembly_directory = sprintf('../assembly_files/order%d_p%d/numRefinements%d/', order, degree, numRefinements);
    
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
    %   Set BCs for the displacement field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Initialize the BC arrays
    BCs_displacement = [];
    BCs_force = [];
    
    
    % Degrees of freedom and BC values at x = -bar_length
    dofIndex = 1;
    bcValue = 0;
    
    elementIndex = 1;
    
    nodeIndex = find_node_index(1, ...
                                [], ...
                                [], ...
                                numNodesPerElement1, ...
                                [], ...
                                [], ...
                                IEN_array(:, elementIndex));
    
    BCs_displacement = [BCs_displacement; ...
                        [GN_array(nodeIndex), dofIndex, bcValue]];
    
    
    % Degrees of freedom and BC values at x = +bar_length
    u_L = 0.1;
    
    dofIndex = 1;
    bcValue = u_L;
    
    elementIndex = numElements1;
    
    nodeIndex = find_node_index(numNodesPerElement1, ...
                                [], ...
                                [], ...
                                numNodesPerElement1, ...
                                [], ...
                                [], ...
                                IEN_array(:, elementIndex));
    
    BCs_displacement = [BCs_displacement; ...
                        [GN_array(nodeIndex), dofIndex, bcValue]];
    
    
    %----------------------------------------------------------------------
    %  Build the ID, LM, and BC arrays
    %----------------------------------------------------------------------
    ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array);
    LM_array1 = build_lm_array(IEN_array, ID_array, GN_array, 0);
    
    [BCU_array1, BCF_array1] = build_bc_array(BCs_displacement, BCs_force, ID_array);
    
    % Number of unknown degrees of freedom
    numUnknownDOFs1 = numDOFs - size(BCU_array1, 1);
    
    % Indices for the unknown coefficients
    index_u1 = (1 : numUnknownDOFs1)';
    index_f1 = ((numUnknownDOFs1 + 1) : numDOFs)';
    
    
    clear BCs_displacement BCs_force ID_array;
    
    
    
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
    
    nodeIndex = find_node_index((1 : (numNodesPerElement1 - 1))', ...
                                [], ...
                                [], ...
                                numNodesPerElement1, ...
                                [], ...
                                [], ...
                                IEN_array(:, elementIndex));
    
    temp = ones(numNodesPerElement1 - 1, 1);
    
    BCs_displacement = [BCs_displacement; ...
                        [GN_array(nodeIndex), dofIndex * temp, bcValue * temp]];
    
    
    % Degrees of freedom and BC values at x = +bar_length
    dofIndex = 1;
    bcValue = 1;
    
    elementIndex = numElements1;
    
    nodeIndex = find_node_index((2 : numNodesPerElement1)', ...
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
                 'LM_array1', 'LM_array2', ...
                 'BCU_array1', 'BCU_array2', ...
                 'BCF_array1', 'BCF_array2', ...
                 'numUnknownDOFs1', 'numUnknownDOFs2', ...
                 'index_u1', 'index_u2', ...
                 'index_f1', 'index_f2', ...
                 'u_L', ...
                 '-v7.3');
end
