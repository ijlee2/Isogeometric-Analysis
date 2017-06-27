%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine builds the ID ("Destination") array from the displacement
%  BC array. A column of the ID array corresponds to a node, and returns
%  the global equation indices for the degrees of freedom (DOFs) on the
%  node.
%  
%  We build the array so that the known displacements (the known B-spline
%  or NURBS coefficients) are placed on the bottom side of the solution
%  vector u. In other words, the unknown coefficients are on the top side.
%  
%  
%  Warning:
%  
%  Use this routine after creating BCs_displacement and BCs_force arrays.
%  These arrays have 3 columns of the form,
%  
%      [global node index, dof index, BC value]
%  
%  Note that we will assume zero force for all DOFs on which displacement
%  or force BC has not been specified. This allows us to specify only the 
%  nonzero forces in the BCs_force array and save memory.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array);
%  
%  where,
%  
%      BCs_displacement is the displacement BC arrays
%      numDOFsPerNode is the number of degree of freedoms that each node
%          has (number of dimensions x number of fields of interest)
%      GN_array is the GN ("Global Node") array
%  
%  
%  Output:
%  
%  1. (numDOFsPerNode) x (numNodes) array
%  
%      ID_array(:, j) returns the global equation indices for the j-th node
%--------------------------------------------------------------------------
function ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array)
    if (~isempty(BCs_displacement))
        %------------------------------------------------------------------
        %  Remove duplicate BC information
        %------------------------------------------------------------------
        BCs_displacement = unique(BCs_displacement, 'rows');
        numBCs_displacement = size(BCs_displacement, 1);
        
        
        % Check for errors
        if (numBCs_displacement ~= size(unique(BCs_displacement(:, [1; 2]), 'rows'), 1))
            fprintf('Warning: Two copies of a node have different DOF values.\n\nNo operation has been performed.\n\n');
            
            ID_array = [];
            
            return;
        end
        
        
        %------------------------------------------------------------------
        %  Remove duplicate BC information due to shared nodes
        %------------------------------------------------------------------
        BCs_displacement = unique([GN_array(BCs_displacement(:, 1)), BCs_displacement(:, 2), BCs_displacement(:, 3)], 'rows');
        numBCs_displacement = size(BCs_displacement, 1);
        
        
        % Check for errors
        if (numBCs_displacement ~= size(unique(BCs_displacement(:, [1; 2]), 'rows'), 1))
            fprintf('Warning: Two nodes that are shared have different DOF values.\n\nNo operation has been performed.\n\n');
            
            ID_array = [];
            
            return;
        end
        
    else
        numBCs_displacement = size(BCs_displacement, 1);
        
    end
    
    
    %----------------------------------------------------------------------
    %  Initialize the ID array
    %----------------------------------------------------------------------
    numNodes = size(GN_array, 1);
    numUniqueNodes = size(unique(GN_array), 1);
    
    ID_array = zeros(numDOFsPerNode, numNodes);
    numDOFs = numDOFsPerNode * numUniqueNodes;
    
    
    %----------------------------------------------------------------------
    %  Assign global equation indices to DOFs that are known
    %----------------------------------------------------------------------
    % Counter for how many displacement DOFs we have encountered
    count = 1;
    
    for i = (numDOFs - numBCs_displacement + 1) : numDOFs
        ID_array(BCs_displacement(count, 2), BCs_displacement(count, 1)) = i;
        
        count = count + 1;
    end
    
    
    %----------------------------------------------------------------------
    %  Assign global equation indices to the remaining DOFs
    %----------------------------------------------------------------------
    % Counter for how many force DOFs we have encountered
    count = 1;
    
    for j = 1 : numNodes
        % If the node is not shared, check whether its DOFs correspond to
        % a displacement BC and already have an equation index assigned
        if (GN_array(j) == j)
            for i = 1 : numDOFsPerNode
                if (ID_array(i, j) == 0)
                    ID_array(i, j) = count;
                    
                    count = count + 1;
                end
            end
            
        % If the node is shared, read the target node's equation indices
        else
            for i = 1 : numDOFsPerNode
                ID_array(i, j) = ID_array(i, GN_array(j));
            end
            
        end
    end
end