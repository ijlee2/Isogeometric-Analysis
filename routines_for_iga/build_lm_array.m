%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine builds the LM ("Location Matrix") array for a B-splines or
%  NURBS patch. A column of the LM array corresponds to an element, and
%  returns the global equation indices for the DOFs in the element.
%  
%  
%  Warning:
%  
%  Use this routine within the constitutive model.
%  
%  We assume that the element matrices will be created such that its rows
%  and columns follow the "usual" DOF order. For example, for 3D quadratic
%  NURBS with one field of interest, we would have,
%  
%     [node 1 - dof 1; ...
%      node 1 - dof 2; ...
%      node 1 - dof 3; ...
%      node 2 - dof 1; ...
%      node 2 - dof 2; ...
%      node 2 - dof 3; ...
%      ...
%      node 27 - dof 1; ...
%      node 27 - dof 2; ...
%      node 27 - dof 3]
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      LM_array = build_lm_array(IEN_array, ID_array, GN_array, numNodesBeforeMe(p));
%  
%  where,
%  
%      IEN_array is the IEN ("Element Nodes") array for the patch
%      ID_array is the ID ("Destination") array
%      GN_array is the GN ("Global Node") array
%      numNodesBeforeMe is the number of nodes that come before the patch
%  
%  
%  Output:
%  
%  1. (numDOFsPerElement) x (numElements) array
%  
%      LM_array(:, e) returns the global equation indices for the e-th
%      element
%--------------------------------------------------------------------------
function LM_array = build_lm_array(IEN_array, ID_array, GN_array, numNodesBeforeMe)
    [numNodesPerElement, numElements] = size(IEN_array);
    numDOFsPerNode = size(ID_array, 1);
    numDOFsPerElement = numDOFsPerNode * numNodesPerElement;
    
    
    %----------------------------------------------------------------------
    %  Initialize the LM array
    %----------------------------------------------------------------------
    LM_array = zeros(numDOFsPerElement, numElements);
    
    
    %----------------------------------------------------------------------
    %  Loop over elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        temp = ID_array(:, GN_array(IEN_array(:, e) + numNodesBeforeMe));
        
        index = 1;
        
        for j = 1 : numNodesPerElement
            for i = 1 : numDOFsPerNode
                LM_array(index, e) = temp(i, j);
                
                index = index + 1;
            end
        end
    end
end