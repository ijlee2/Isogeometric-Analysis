%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  For a B-spline or NURBS patch, this routine finds the global indices of
%  the elements that match the local element indices along each direction.
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
%      elementIndex = find_element_index(elementLocalIndex1, elementLocalIndex2, elementLocalIndex3, numElements1, numElements2, numElements3);  (for 3D)
%      elementIndex = find_element_index(elementLocalIndex1, elementLocalIndex2, [], numElements1, numElements2, []);  (for 2D)
%      elementIndex = find_element_index(elementLocalIndex1, [], [], numElements1, [], []);  (for 1D)
%  
%  where,
%  
%      elementLocalIndex1, elementLocalIndex2, elementLocalIndex3 are
%          arrays of local element indices along directions 1, 2, and 3
%          (column vectors)
%      numElements1, numElements2, numElements3 are the numbers of elements
%          along directions 1, 2, and 3
%  
%  
%  Output:
%  
%  1. (numElementsSpecified) x 1 column vector
%  
%      The column contains the global indices of the elements
%--------------------------------------------------------------------------
function elementIndex = find_element_index(elementLocalIndex1, elementLocalIndex2, elementLocalIndex3, numElements1, numElements2, numElements3)
    numElementsSpecified1 = size(elementLocalIndex1, 1);
    numElementsSpecified2 = size(elementLocalIndex2, 1);
    numElementsSpecified3 = size(elementLocalIndex3, 1);
    
    
    % Check if the element is in 1D
    if (numElementsSpecified2 == 0 && numElementsSpecified3 == 0)
        % Set the numbers of elements
        numElements2 = 1;
        numElements3 = 1;
        
        % Set the local element index arrays to be 1 x 1 vectors of ones
        elementLocalIndex2 = 1;
        elementLocalIndex3 = 1;
        numElementsSpecified2 = 1;
        numElementsSpecified3 = 1;
        
    % Check if the element is in 2D
    elseif (numElementsSpecified3 == 0)
        % Set the number of elements
        numElements3 = 1;
        
        % Set the local element index array to be a 1 x 1 vector of ones
        elementLocalIndex3 = 1;
        numElementsSpecified3 = 1;
        
    end
    
    
    % Check for errors
    if (min(elementLocalIndex1) < 1 || max(elementLocalIndex1) > numElements1)
        fprintf('There are %d elements along direction 1, so the local element index array\nfor direction 1 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numElements1, numElements1);
        
        elementIndex = [];
        
        return;
    end
    
    if (min(elementLocalIndex2) < 1 || max(elementLocalIndex2) > numElements2)
        fprintf('There are %d elements along direction 2, so the local element index array\nfor direction 2 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numElements2, numElements2);
        
        elementIndex = [];
        
        return;
    end
    
    if (min(elementLocalIndex3) < 1 || max(elementLocalIndex3) > numElements3)
        fprintf('There are %d elements along direction 3, so the local element index array\nfor direction 3 can only contain numbers between 1 and %d.\n\nNo operation has been performed.\n\n', numElements3, numElements3);
        
        elementIndex = [];
        
        return;
    end
    
    
    %----------------------------------------------------------------------
    %  Find the global element indices
    %----------------------------------------------------------------------
%   numElementsSpecified = numElementsSpecified1 * numElementsSpecified2 * numElementsSpecified3;
    
    
    % Local node element contributions from direction 1, 2, and 3
    temp3 = kron((numElements1 * numElements2) * (elementLocalIndex3 - 1), ones(numElementsSpecified1 * numElementsSpecified2, 1));
    
    temp2 = repmat(kron(numElements1 * (elementLocalIndex2 - 1), ones(numElementsSpecified1, 1)), numElementsSpecified3, 1);
    
    temp1 = repmat(elementLocalIndex1, numElementsSpecified2 * numElementsSpecified3, 1);
    
    elementIndex = temp1 + temp2 + temp3;
end