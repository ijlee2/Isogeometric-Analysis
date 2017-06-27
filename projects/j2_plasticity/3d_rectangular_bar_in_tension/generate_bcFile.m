function generate_bcFile(numRefinements)
    % Set the path to the assembly files directory
    path_to_assembly_directory = sprintf('../assembly_files/numRefinements%d/', numRefinements);
    
    % File paths
    file_assembly_global = strcat(path_to_assembly_directory, 'file_assembly_global');
    file_assembly_patch  = strcat(path_to_assembly_directory, 'file_assembly_patch');
    
    % Load the global assembly file
    load(file_assembly_global, 'GN_array');
    
    % Load patch 1
    load(sprintf('%s%d', file_assembly_patch, 1), ...
         'numNodesPerElement1', ...
         'numNodesPerElement2', ...
         'numNodesPerElement3', ...
         'numElements1', ...
         'numElements2', ...
         'numElements3', ...
         'IEN_array');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Set load step 1
    %   
    %   Here, we fix the top grip and displace the bottom grip in the
    %   x-direction by 10 cm.
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Set the current load step index
    loadStep = 1;
    
    % Set the number of boundaries
    numBoundaries = 2;
    
    % Set the default number of time steps
    numTimeSteps = 1000;
    
    % Save load step 1
    save(sprintf('%sfile_loadstep%d', path_to_assembly_directory, loadStep), ...
                 'numBoundaries', ...
                 'numTimeSteps', ...
                 '-v7.3');
    
    
    %----------------------------------------------------------------------
    %  Set the nodes on boundary 1 (top grip)
    %----------------------------------------------------------------------
    % Set the current boundary index
    boundary = 1;
    
    % Initialize the nodes on the boundary
    nodesOnBoundary = [];
    
    
    % Find the element index of the grip along direction 1
    switch numRefinements
        case 0
            numElements1_grip = 2;
            elementIndex1_grip = (1 : numElements1_grip)';
            
        case 1
            numElements1_grip = 4;
            elementIndex1_grip = (1 : numElements1_grip)';
        
        case 2
            numElements1_grip = 4;
            elementIndex1_grip = (1 : numElements1_grip)';
        
    end
    
    
    % y-z surface, x = 0
    elementIndex = find_element_index(1, (1 : numElements2)', (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index(1, (1 : numNodesPerElement2)', (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-y surface, z = 0
    elementIndex = find_element_index(elementIndex1_grip, (1 : numElements2)', 1, numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', (1 : numNodesPerElement2)', 1, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-y surface, z = bar_thickness
    elementIndex = find_element_index(elementIndex1_grip, (1 : numElements2)', numElements3, numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', (1 : numNodesPerElement2)', numNodesPerElement3, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-z surface, y = 0
    elementIndex = find_element_index(elementIndex1_grip, 1, (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', 1, (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-z surface, y = grip_width
    elementIndex = find_element_index(elementIndex1_grip, numElements2, (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', numNodesPerElement2, (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    %----------------------------------------------------------------------
    %  Set the loads on boundary 1 (top grip)
    %  
    %  Load type = 1, for displacement only
    %            = 2, for rotation only
    %            = 3, for both displacement and rotation
    %----------------------------------------------------------------------
    % Set the load type
    loadType = 1;
    
    % Remove the duplicate nodes
    nodesOnBoundary = unique(nodesOnBoundary);
    
    % Set the degrees of freedom that are affected
    dofIndex = [1; 2; 3];
    
    
    % Specify the total displacement vector
    displacementTotal = [0; 0; 0];
    
    % Specify the total rotation
    rotationTotal = 0;
    
    % Specify the unit vector about which rotation occurs
    rotationAxis = [1; 0; 0];
    
    
    %----------------------------------------------------------------------
    %  Save load step 1 for boundary 1
    %----------------------------------------------------------------------
    save(sprintf('%sfile_loadstep%d_boundary%d', path_to_assembly_directory, loadStep, boundary), ...
                 'loadType', ...
                 'nodesOnBoundary', ...
                 'dofIndex', ...
                 'displacementTotal', ...
                 'rotationTotal', ...
                 'rotationAxis', ...
                 '-v7.3');
    
    
    %----------------------------------------------------------------------
    %  Set the nodes on boundary 2 (bottom grip)
    %----------------------------------------------------------------------
    % Set the current boundary index
    boundary = 2;
    
    % Initialize the nodes on the boundary
    nodesOnBoundary = [];
    
    
    % Find the element index of the grip along direction 1
    switch numRefinements
        case 0
            numElements1_grip = 2;
            elementIndex1_grip = ((numElements1 - numElements1_grip + 1) : numElements1)';
            
        case 1
            numElements1_grip = 4;
            elementIndex1_grip = ((numElements1 - numElements1_grip + 1) : numElements1)';
        
        case 2
            numElements1_grip = 4;
            elementIndex1_grip = ((numElements1 - numElements1_grip + 1) : numElements1)';
            
    end
    
    
    % y-z surface, x = bar_totalLength
    elementIndex = find_element_index(numElements1, (1 : numElements2)', (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index(numNodesPerElement1, (1 : numNodesPerElement2)', (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-y surface, z = 0
    elementIndex = find_element_index(elementIndex1_grip, (1 : numElements2)', 1, numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', (1 : numNodesPerElement2)', 1, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-y surface, z = bar_thickness
    elementIndex = find_element_index(elementIndex1_grip, (1 : numElements2)', numElements3, numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', (1 : numNodesPerElement2)', numNodesPerElement3, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-z surface, y = 0
    elementIndex = find_element_index(elementIndex1_grip, 1, (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', 1, (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    % x-z surface, y = grip_width
    elementIndex = find_element_index(elementIndex1_grip, numElements2, (1 : numElements3)', numElements1, numElements2, numElements3);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', numNodesPerElement2, (1 : numNodesPerElement3)', numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, elementIndex));
    
    nodesOnBoundary = [nodesOnBoundary; GN_array(nodeIndex)];
    
    
    %----------------------------------------------------------------------
    %  Set the loads on boundary 2 (top grip)
    %  
    %  Load type = 1, for displacement only
    %            = 2, for rotation only
    %            = 3, for both displacement and rotation
    %----------------------------------------------------------------------
    % Set the load type
    loadType = 1;
    
    % Remove the duplicate nodes
    nodesOnBoundary = unique(nodesOnBoundary);
    
    % Set the degrees of freedom that are affected
    dofIndex = [1; 2; 3];
    
    
    % Specify the total displacement vector
    displacementTotal = [0.05; 0; 0];
    
    % Specify the total rotation
    rotationTotal = 0;
    
    % Specify the unit vector about which rotation occurs
    rotationAxis = [1; 0; 0];
    
    
    %----------------------------------------------------------------------
    %  Save load step 2 for boundary 2
    %----------------------------------------------------------------------
    save(sprintf('%sfile_loadstep%d_boundary%d', path_to_assembly_directory, loadStep, boundary), ...
                 'loadType', ...
                 'nodesOnBoundary', ...
                 'dofIndex', ...
                 'displacementTotal', ...
                 'rotationTotal', ...
                 'rotationAxis', ...
                 '-v7.3');
end
