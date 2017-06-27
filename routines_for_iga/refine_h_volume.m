%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine does a h-refinement (knot refinement) of the given knot
%  vectors in 3D.
%  
%  
%  Warning:
%  
%  The knot vectors are assumed to be open, and the knots for insertion are
%  assumed to be in the interior of the knot vector. The nodes are assumed
%  to be ordered along the direction of 1, then along that of 2, then along
%  that of 3.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [knots1_new, knots2_new, knots3_new, nodes_new] = refine_h_volume(knots1, knots2, knots3, nodes, p1, p2, p3, knotsForInsertion1, knotsForInsertion2, knotsForInsertion3);
%  
%  where,
%  
%      knots1, knots2, knots3 are the knot vectors (column vectors)
%      nodes is an array of nodes
%      p1, p2, p3 are the degree of the B-splines
%      knotsForInsertion1, knotsForInsertion2, knotsForInsertion3 are the
%          knot vectors that we seek to insert (column vectors)
%  
%  
%  Outputs:
%  
%  1. New knot vectors (same degrees)
%  
%  2. New nodes array
%--------------------------------------------------------------------------
function [knots1_new, knots2_new, knots3_new, nodes_new] = refine_h_volume(knots1, knots2, knots3, nodes, p1, p2, p3, knotsForInsertion1, knotsForInsertion2, knotsForInsertion3)
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    constant_p3p1 = p3 + 1;
    
    numKnots1 = size(knots1, 1);
    numKnots2 = size(knots2, 1);
    numKnots3 = size(knots3, 1);
    numNodes1 = numKnots1 - constant_p1p1;
    numNodes2 = numKnots2 - constant_p2p1;
    numNodes3 = numKnots3 - constant_p3p1;
    numNodes  = numNodes1 * numNodes2 * numNodes3;
    numDimensions = size(nodes, 2);
    
    numKnotInsertions1 = size(knotsForInsertion1, 1);
    numKnotInsertions2 = size(knotsForInsertion2, 1);
    numKnotInsertions3 = size(knotsForInsertion3, 1);
    numKnots1_new = numKnots1 + numKnotInsertions1;
    numKnots2_new = numKnots2 + numKnotInsertions2;
    numKnots3_new = numKnots3 + numKnotInsertions3;
    numNodes1_new = numNodes1 + numKnotInsertions1;
    numNodes2_new = numNodes2 + numKnotInsertions2;
    numNodes3_new = numNodes3 + numKnotInsertions3;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots1_new = knots1;
        knots2_new = knots2;
        knots3_new = knots3;
        nodes_new = nodes;
        
        return;
        
    elseif (numKnotInsertions1 == 0 && numKnotInsertions2 == 0 && numKnotInsertions3 == 0)
        fprintf('Warning: Please specify at least one knot for insertion.\n\nNo operation has been performed.\n\n');
        
        knots1_new = knots1;
        knots2_new = knots2;
        knots3_new = knots3;
        nodes_new = nodes;
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Knot refinement for direction 1
    %----------------------------------------------------------------------
    if (numKnotInsertions1 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "plane"
        %------------------------------------------------------------------
        temp2 = (0 : (numNodes2 - 1))';
        temp3 = (0 : (numNodes3 - 1))';
        ones2 = ones(numNodes2, 1);
        
        shifts_for_nodes     = repmat(numNodes1     * temp2, numNodes3, 1) + kron((numNodes1     * numNodes2) * temp3, ones2);
        shifts_for_nodes_new = repmat(numNodes1_new * temp2, numNodes3, 1) + kron((numNodes1_new * numNodes2) * temp3, ones2);
        
        clear temp2 temp3 ones2;
        
        
        %------------------------------------------------------------------
        %  Initialize the new knot vector and nodes array
        %------------------------------------------------------------------
        knots1_new = [knots1; zeros(numKnotInsertions1, 1)];
        
        clear knots1;
        
        nodes_new = zeros(numNodes1_new * numNodes2 * numNodes3, numDimensions);
        
        for i = 1 : numNodes1
            % For the last node along direction 1, we make extra copies
            if (i < numNodes1)
                nodes_new(i + shifts_for_nodes_new, :) = nodes(i + shifts_for_nodes, :);
            else
                for j = 0 : numKnotInsertions1
                    nodes_new(i + j + shifts_for_nodes_new, :) = nodes(i + shifts_for_nodes, :);
                end
            end
        end
        
        clear nodes shifts_for_nodes;
        
        
        %------------------------------------------------------------------
        %  Loop over the knots for insertion
        %------------------------------------------------------------------
        % Search where the knot xi belongs at this index
        knotSearch_begin = constant_p1p1;
        
        for b = 1 : numKnotInsertions1
            %--------------------------------------------------------------
            %  Knot insertion
            %--------------------------------------------------------------
            % Knot that we will insert
            xi = knotsForInsertion1(b);
            
            % Find where to insert the knot
            for k = knotSearch_begin : (numKnots1_new - 1)
                if (knots1_new(k) <= xi && xi < knots1_new(k + 1))
                    index_knot = k + 1;
                    
                    break;
                end
            end
            
            % Shift the "end" nodes
            for i = (numNodes1 + b - 1) : -1 : index_knot
                nodes_new(i + shifts_for_nodes_new, :) = nodes_new(i - 1 + shifts_for_nodes_new, :);
            end
            
            % Update the "middle" nodes
            for i = (index_knot - 1) : -1 : (index_knot - p1)
                alpha = (xi - knots1_new(i)) / (knots1_new(i + p1) - knots1_new(i));
                
                nodes_new(i + shifts_for_nodes_new, :) = alpha  * nodes_new(i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(i - 1 + shifts_for_nodes_new, :);
            end
            
            % Shift the "end" knots
            for i = (numKnots1 + b) : -1 : (index_knot + 1)
                knots1_new(i) = knots1_new(i - 1);
            end
            
            
            %--------------------------------------------------------------
            %  Insert knot xi to the new knot vector
            %--------------------------------------------------------------
            knots1_new(index_knot) = xi;
            
            
            %--------------------------------------------------------------
            %  Update index for the next iteration
            %--------------------------------------------------------------
            knotSearch_begin = index_knot;
        end
        
        nodes = nodes_new;
        
        clear nodes_new;
        
    % Default action for no refinement
    else
        knots1_new = knots1;
        
        clear knots1;
    end
    
    
    %----------------------------------------------------------------------
    %  Knot refinement for direction 2
    %----------------------------------------------------------------------
    if (numKnotInsertions2 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "plane"
        %  
        %  Note that numNodes1 = numNodes1_new now. We add an extra shift
        %  of -(numNodes1_new - 1) so that we can sequentially traverse
        %  along direction 2 like we did along direction 1.
        %------------------------------------------------------------------
        temp1 = (0 : numNodes1_new - 1)';
        temp3 = (0 : numNodes3 - 1)';
        ones1 = ones(numNodes1_new, 1);
        
        shifts_for_nodes     = repmat(temp1, numNodes3, 1) + kron((numNodes1_new * numNodes2    ) * temp3, ones1) - (numNodes1_new - 1);
        shifts_for_nodes_new = repmat(temp1, numNodes3, 1) + kron((numNodes1_new * numNodes2_new) * temp3, ones1) - (numNodes1_new - 1);
        
        clear temp1 temp3 ones1;
        
        
        %------------------------------------------------------------------
        %  Initialize the new knot vector and nodes array
        %------------------------------------------------------------------
        knots2_new = [knots2; zeros(numKnotInsertions2, 1)];
        
        clear knots2;
        
        nodes_new = zeros(numNodes1 * numNodes2_new * numNodes3, numDimensions);
        
        for i = 1 : numNodes2
            % For the last node along direction 1, we make extra copies
            if (i < numNodes2)
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = nodes(numNodes1_new * i + shifts_for_nodes, :);
            else
                for j = 0 : numKnotInsertions2
                    nodes_new(numNodes1_new * (i + j) + shifts_for_nodes_new, :) = nodes(numNodes1_new * i + shifts_for_nodes, :);
                end
            end
        end
        
        clear nodes shifts_for_nodes;
        
        
        %------------------------------------------------------------------
        %  Loop over the knots for insertion
        %------------------------------------------------------------------
        % Search where the knot xi belongs at this index
        knotSearch_begin = constant_p2p1;
        
        for b = 1 : numKnotInsertions2
            %--------------------------------------------------------------
            %  Knot insertion
            %--------------------------------------------------------------
            % Knot that we will insert
            xi = knotsForInsertion2(b);
            
            % Find where to insert the knot
            for k = knotSearch_begin : (numKnots2_new - 1)
                if (knots2_new(k) <= xi && xi < knots2_new(k + 1))
                    index_knot = k + 1;
                    
                    break;
                end
            end
            
            % Shift the "end" nodes
            for i = (numNodes2 + b - 1) : -1 : index_knot
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = nodes_new(numNodes1_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Update the "middle" nodes
            for i = (index_knot - 1) : -1 : (index_knot - p2)
                alpha = (xi - knots2_new(i)) / (knots2_new(i + p2) - knots2_new(i));
                
                nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = alpha * nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(numNodes1_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Shift the "end" knots
            for i = (numKnots2 + b) : -1 : (index_knot + 1)
                knots2_new(i) = knots2_new(i - 1);
            end
            
            
            %--------------------------------------------------------------
            %  Insert knot xi to the new knot vector
            %--------------------------------------------------------------
            knots2_new(index_knot) = xi;
            
            
            %--------------------------------------------------------------
            %  Update index for the next iteration
            %--------------------------------------------------------------
            knotSearch_begin = index_knot;
        end
        
        nodes = nodes_new;
        
        clear nodes_new;
        
    % Default action for no refinement
    else
        knots2_new = knots2;
        
        clear knots2;
    end
    
    
    %----------------------------------------------------------------------
    %  Knot refinement for direction 3
    %----------------------------------------------------------------------
    if (numKnotInsertions3 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "plane"
        %  
        %  We have numNodes1 = numNodes1_new and numNodes2 = numNodes2_new.
        %  We add an extra shift of -(numNodes1_new * numNodes2_new - 1)
        %  so that we can sequentially traverse along direction 3 like we
        %  did along direction 1.
        %------------------------------------------------------------------
        % Number of nodes on the 1-2 plane
        numNodes12_new = numNodes1_new * numNodes2_new;
        
        shifts_for_nodes     = (0 : numNodes12_new - 1)' - (numNodes12_new - 1);
        shifts_for_nodes_new = shifts_for_nodes;
        
        
        %------------------------------------------------------------------
        %  Initialize the new knot vector and nodes array
        %------------------------------------------------------------------
        knots3_new = [knots3; zeros(numKnotInsertions3, 1)];
        
        clear knots3;
        
        nodes_new = zeros(numNodes12_new * numNodes3_new, numDimensions);
        
        for i = 1 : numNodes3
            % For the last node along direction 3, we make extra copies
            if (i < numNodes3)
                nodes_new(numNodes12_new * i + shifts_for_nodes_new, :) = nodes(numNodes12_new * i + shifts_for_nodes, :);
            else
                for j = 0 : numKnotInsertions3
                    nodes_new(numNodes12_new * (i + j) + shifts_for_nodes_new, :) = nodes(numNodes12_new * i + shifts_for_nodes, :);
                end
            end
        end
        
        clear nodes shifts_for_nodes;
        
        
        %------------------------------------------------------------------
        %  Loop over the knots for insertion
        %------------------------------------------------------------------
        % Search where the knot xi belongs at this index
        knotSearch_begin = constant_p3p1;
        
        for b = 1 : numKnotInsertions3
            %--------------------------------------------------------------
            %  Knot insertion
            %--------------------------------------------------------------
            xi = knotsForInsertion3(b);
            
            % Find where to insert the knot
            for k = knotSearch_begin : (numKnots3_new - 1)
                if (knots3_new(k) <= xi && xi < knots3_new(k + 1))
                    index_knot = k + 1;
                    
                    break;
                end
            end
            
            % Shift the "end" nodes
            for i = (numNodes3 + b - 1) : -1 : index_knot
                nodes_new(numNodes12_new * i + shifts_for_nodes_new, :) = nodes_new(numNodes12_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Update the "middle" nodes
            for i = (index_knot - 1) : -1 : (index_knot - p3)
                alpha = (xi - knots3_new(i)) / (knots3_new(i + p3) - knots3_new(i));
                
                nodes_new(numNodes12_new * i + shifts_for_nodes_new, :) = alpha * nodes_new(numNodes12_new * i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(numNodes12_new * (i - 1) + shifts_for_nodes_new, :);
            end
            
            % Shift the "end" knots
            for i = (numKnots3 + b) : -1 : (index_knot + 1)
                knots3_new(i) = knots3_new(i - 1);
            end
            
            
            %--------------------------------------------------------------
            %  Insert knot xi to the new knot vector
            %--------------------------------------------------------------
            knots3_new(index_knot) = xi;
            
            
            %--------------------------------------------------------------
            %  Update index for the next iteration
            %--------------------------------------------------------------
            knotSearch_begin = index_knot;
        end
        
    % Default action for no refinement
    else
        knots3_new = knots3;
        nodes_new = nodes;
        
        clear knots3 nodes;
    end
end