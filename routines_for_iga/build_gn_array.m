%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine builds the GN ("Global Node") array for B-spline or NURBS
%  patches.
%  
%  
%  Warning:
%  
%  SN_array has 4 columns of the form,
%  
%      [my patch index, my node index, target patch index, target node index]
%  
%  When linking a node (my node) to another node (target node), we respect
%  the order in which the patches have been created and the "usual" way
%  by which we march along the nodes. Thus, we always link a node to one
%  (and only one) node that comes beforehand.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      GN_array = build_gn_array(SN_array, numNodes);
%  
%  where,
%  
%      SN_array is the shared nodes array
%      numNodesOnPatch is an array of the number of nodes on each patch
%          (column vector)
%  
%  
%  Output:
%  
%  1. (numNodes) x 1 vector
%  
%      GN_array(i) returns the index of the node that that i-th node
%      coincides with
%--------------------------------------------------------------------------
function GN_array = build_gn_array(SN_array, numNodesOnPatch)
    % Number of patches
    numPatches = size(numNodesOnPatch, 1);
    
    % Number of nodes created before the current patch
    numNodesBeforeMe = zeros(numPatches, 1);
    
    % Total number of nodes
    numNodes = numNodesOnPatch(1);
    
    for p = 2 : numPatches
        numNodesBeforeMe(p) = numNodesBeforeMe(p - 1) + numNodesOnPatch(p - 1);
        numNodes = numNodes + numNodesOnPatch(p);
    end
    
    
    %----------------------------------------------------------------------
    %  Initialize the GN array
    %----------------------------------------------------------------------
    % By default, we assume that none of the nodes coincide with one another
    GN_array = (1 : numNodes)';
    
    
    %----------------------------------------------------------------------
    %  Loop over shared nodes
    %----------------------------------------------------------------------
    [SN_array, index_original] = unique(SN_array, 'rows');
    
    % Number of nodes on the patch that coincide with other nodes
    numSharedNodes = size(SN_array, 1);
    
    
    for i = 1 : numSharedNodes
        % Check for errors
        if (SN_array(i, 3) > SN_array(i, 1))
            fprintf('Warning: Check row %d of shared nodes array.\n\nTarget patch index cannot exceed the current patch index.\n\nNo operation has been performed.\n\n', index_original(i));
            
            GN_array = [];
            
            return;
            
        elseif (SN_array(i, 3) == SN_array(i, 1) && SN_array(i, 4) > SN_array(i, 2))
            fprintf('Warning: Check row %d of shared nodes array.\n\nTarget node index cannot exceed the current node index if they are from the same patch.\n\nNo operation has been performed.\n\n', index_original(i));
            
            GN_array = [];
            
            return;
            
        elseif (SN_array(i, 2) > numNodesOnPatch(SN_array(i, 1)))
            fprintf('Warning: Check row %d of shared nodes array.\n\nCurrent node index cannot exceed the number of nodes on the current patch.\n\nNo operation has been performed.\n\n', index_original(i));
            
            GN_array = [];
            
            return;
            
        elseif (SN_array(i, 4) > numNodesOnPatch(SN_array(i, 3)))
            fprintf('Warning: Check row %d of shared nodes array.\n\nTarget node index cannot exceed the number of nodes on the target patch.\n\nNo operation has been performed.\n\n', index_original(i));
            
            GN_array = [];
            
            return;
            
        elseif (i < numSharedNodes && (SN_array(i, 1) == SN_array(i + 1, 1) && SN_array(i, 2) == SN_array(i + 1, 2)))
            fprintf('Warning: Check rows %d and %d of shared nodes array.\n\nPlease assign the current node to only one target node.\n\nNo operation has been performed.\n\n', index_original(i), index_original(i + 1));
            
            GN_array = [];
            
            return;
            
        end
        
        
        % My node index
        myNodeIndex = numNodesBeforeMe(SN_array(i, 1)) + SN_array(i, 2);
        
        % Target node index
        targetNodeIndex = numNodesBeforeMe(SN_array(i, 3)) + SN_array(i, 4);
        
        % Find the "root" of the target node
        while (GN_array(targetNodeIndex) ~= targetNodeIndex)
            targetNodeIndex = GN_array(targetNodeIndex);
        end
        
        % Set the global node index to the target global node index
        GN_array(myNodeIndex) = GN_array(targetNodeIndex);
    end
end