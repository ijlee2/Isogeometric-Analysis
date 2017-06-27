%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine builds the IEN ("Element Nodes") array for a B-splines or
%  NURBS patch. A column of the IEN array corresponds to an element,
%  and returns the global indices of the nodes that make up the element.
%  
%  
%  Warning:
%  
%  Use this routine after calling build_bezier_extraction to get the node
%  index shifts arrays. The knot vectors are assumed to be open.
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,%  
%      IEN_array = build_ien_array(knots1, knots2, knots3, p1, p2, p3, nodeIndexShifts1, nodeIndexShifts2, nodeIndexShifts3);  (for 3D)
%      IEN_array = build_ien_array(knots1, knots2, [], p1, p2, [], nodeIndexShifts1, nodeIndexShifts2, []);                    (for 2D)
%      IEN_array = build_ien_array(knots1, [], [], p1, [], [], nodeIndexShifts1, [], []);                                      (for 1D)
%  
%  where,
%  
%      knots1, knots2, knots3 are the knot vectors (column vectors)
%      p1, p2, p3 are the degrees of the B-splines
%      nodeIndexShifts1, nodeIndexShifts2, nodeIndexShifts3 are the node
%          index shift arrays
%  
%  
%  Output:
%  
%  1. (numNodesPerElement) x (numElements) array
%  
%      IEN_array(:, e) returns the global node indices for the e-th element
%--------------------------------------------------------------------------
function IEN_array = build_ien_array(knots1, knots2, knots3, p1, p2, p3, nodeIndexShifts1, nodeIndexShifts2, nodeIndexShifts3)
    numKnots1 = size(knots1, 1);
    numKnots2 = size(knots2, 1);
    numKnots3 = size(knots3, 1);
    
    % Check if the patch is in 1D
    if (numKnots2 == 0 && numKnots3 == 0)
        % Set the degrees to be 0
        p2 = 0;
        p3 = 0;
        
        % Set the node index shifts arrays to be 1 x 1 zero vectors
        nodeIndexShifts2 = 0;
        nodeIndexShifts3 = 0;
    
    % Check if the patch is in 2D
    elseif (numKnots3 == 0)
        % Set the degree to be 0
        p3 = 0;
        
        % Set the node index shifts array to be a 1 x 1 zero vector
        nodeIndexShifts3 = 0;
        
    end
    
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    constant_p3p1 = p3 + 1;
    
    numNodes1 = numKnots1 - constant_p1p1;
    numNodes2 = numKnots2 - constant_p2p1;
%   numNodes3 = numKnots3 - constant_p3p1;
    
    numElements1 = size(nodeIndexShifts1, 1);
    numElements2 = size(nodeIndexShifts2, 1);
    numElements3 = size(nodeIndexShifts3, 1);
    
    
    
    %----------------------------------------------------------------------
    %  Initialize the IEN array
    %----------------------------------------------------------------------
    numNodesPerElement = constant_p1p1 * constant_p2p1 * constant_p3p1;
    numElements = numElements1 * numElements2 * numElements3;
    
    IEN_array = zeros(numNodesPerElement, numElements);
    
    
    %----------------------------------------------------------------------
    %  Loop over elements
    %----------------------------------------------------------------------
    % Counter for the element
    e = 1;
    
    for e3 = 1 : numElements3
        % Node index contribution from direction 3
        temp3 = kron((numNodes1 * numNodes2) * (((e3 - 1) : (e3 + p3 - 1))' + nodeIndexShifts3(e3)), ones(constant_p1p1 * constant_p2p1, 1));
        
        for e2 = 1 : numElements2
            % Node index contribution from direction 2
            temp2 = repmat(kron(numNodes1 * (((e2 - 1) : (e2 + p2 - 1))' + nodeIndexShifts2(e2)), ones(constant_p1p1, 1)), constant_p3p1, 1);
            
            for e1 = 1 : numElements1
                % Node index contribution from direction 1
                temp1 = repmat((e1 : (e1 + p1))' + nodeIndexShifts1(e1), constant_p2p1 * constant_p3p1, 1);
                
                % Store the indices of the nodes that make up the element
                IEN_array(:, e) = temp1 + temp2 + temp3;
                
                % Increment the counter
                e = e + 1;
            end
        end
    end
end