%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  For a B-spline or NURBS element, this routine finds the global indices
%  of the nodes that match the local node indices along each direction.
%  
%  
%  Warning:
%  
%  Use this routine after calling build_ien_array routine. Note that we
%  can pass as many column(s) of the IEN array as we like. This is useful
%  if we want to prescribe a BC over multiple elements "on the same side."
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,
%  
%      nodeIndex = find_node_index(nodeLocalIndex1, nodeLocalIndex2, nodeLocalIndex3, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array(:, e));  (for 3D)
%      nodeIndex = find_node_index(nodeLocalIndex1, nodeLocalIndex2, [], numNodesPerElement1, numNodesPerElement2, [], IEN_array(:, e));  (for 2D)
%      nodeIndex = find_node_index(nodeLocalIndex1, [], [] , numNodesPerElement1, [], [], IEN_array(:, e));  (for 1D)
%  
%  where,
%  
%      nodeLocalIndex1, nodeLocalIndex2, nodeLocalIndex3 are arrays of
%          local node indices along directions 1, 2, and 3 (column vectors)
%      numNodesPerElement1, numNodesPerElement2, numNodesPerElement3 are
%          the number of nodes on the element along directions 1, 2, and 3
%          (they equal to p1 + 1, p2 + 1, and p3 + 1)
%      IEN_array is the IEN ("Element Nodes") array
%  
%  
%  Output:
%  
%  1. (numNodesSpecified * numElementsSpecified) x 1 column vector
%  
%      The column contains the global indices of the nodes on the elements
%      that were specified
%--------------------------------------------------------------------------
function nodeIndex = find_node_index(nodeLocalIndex1, nodeLocalIndex2, nodeLocalIndex3, numNodesPerElement1, numNodesPerElement2, numNodesPerElement3, IEN_array_e)
    numNodesSpecified1 = size(nodeLocalIndex1, 1);
    numNodesSpecified2 = size(nodeLocalIndex2, 1);
    numNodesSpecified3 = size(nodeLocalIndex3, 1);
    
    
    % Check if the element is in 1D
    if (numNodesSpecified2 == 0 && numNodesSpecified3 == 0)
        % Set the degrees
        numNodesPerElement2 = 1;
        numNodesPerElement3 = 1;
        
        % Set the local node index arrays to be 1 x 1 vectors of ones
        nodeLocalIndex2 = 1;
        nodeLocalIndex3 = 1;
        numNodesSpecified2 = 1;
        numNodesSpecified3 = 1;
        
    % Check if the element is in 2D
    elseif (numNodesSpecified3 == 0)
        % Set the degree
        numNodesPerElement3 = 1;
        
        % Set the local node index array to be a 1 x 1 vector of ones
        nodeLocalIndex3 = 1;
        numNodesSpecified3 = 1;
        
    end
    
    
    % Check for errors
    if (min(nodeLocalIndex1) < 1 || max(nodeLocalIndex1) > numNodesPerElement1)
        fprintf('There are %d nodes along direction 1, so the local node index array\nfor direction 1 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numNodesPerElement1, numNodesPerElement1);
        
        nodeIndex = [];
        
        return;
    end
    
    if (min(nodeLocalIndex2) < 1 || max(nodeLocalIndex2) > numNodesPerElement2)
        fprintf('There are %d nodes along direction 2, so the local node index array\nfor direction 2 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numNodesPerElement2, numNodesPerElement2);
        
        nodeIndex = [];
        
        return;
    end
    
    if (min(nodeLocalIndex3) < 1 || max(nodeLocalIndex3) > numNodesPerElement3)
        fprintf('There are %d nodes along direction 3, so the local node index array\nfor direction 3 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numNodesPerElement3, numNodesPerElement3);
        
        nodeIndex = [];
        
        return;
    end
    
    
    %----------------------------------------------------------------------
    %  Find the global node indices
    %----------------------------------------------------------------------
    numNodesSpecified = numNodesSpecified1 * numNodesSpecified2 * numNodesSpecified3;
    numElementsSpecified = size(IEN_array_e, 2);
    
    
    % Local node index contributions from direction 1, 2, and 3
    temp3 = kron((numNodesPerElement1 * numNodesPerElement2) * (nodeLocalIndex3 - 1), ones(numNodesSpecified1 * numNodesSpecified2, 1));
    
    temp2 = repmat(kron(numNodesPerElement1 * (nodeLocalIndex2 - 1), ones(numNodesSpecified1, 1)), numNodesSpecified3, 1);
    
    temp1 = repmat(nodeLocalIndex1, numNodesSpecified2 * numNodesSpecified3, 1);
    
    % Use the IEN array to find the global node indices
    if (numElementsSpecified == 1)
        nodeIndex = IEN_array_e(temp1 + temp2 + temp3, :);
    else
        nodeIndex = reshape(IEN_array_e(temp1 + temp2 + temp3, :), numNodesSpecified * numElementsSpecified, 1);
    end
end