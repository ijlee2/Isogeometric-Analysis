%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine does a p-refinement (degree elevation) of the given knot
%  knot vector in 2D.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open. The nodes are assumed to be
%  ordered along the direction of 1, then along that of 2.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [knots1_new, knots2_new, nodes_new, p1_new, p2_new] = refine_p_surface(knots1, knots2, nodes, p1, p2, q1, q2);
%  
%  where,
%  
%      knots1, knots2 are the knot vectors (column vectors)
%      nodes is an array of nodes
%      p1, p2 are the degree of the B-splines
%      q1, q2 are the degree increments
%  
%  
%  Outputs:
%  
%  1. New knot vectors
%  
%  2. New nodes array
%  
%  3. New degrees
%--------------------------------------------------------------------------
function [knots1_new, knots2_new, nodes_new, p1_new, p2_new] = refine_p_surface(knots1, knots2, nodes, p1, p2, q1, q2)
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    constant_p1pq1 = p1 + q1;
    constant_p2pq2 = p2 + q2;
    constant_p1pq1_p1 = constant_p1pq1 + 1;
    constant_p2pq2_p1 = constant_p2pq2 + 1;
    
    numKnots1 = size(knots1, 1);
    numKnots2 = size(knots2, 1);
    numNodes1 = numKnots1 - constant_p1p1;
    numNodes2 = numKnots2 - constant_p2p1;
    numNodes  = numNodes1 * numNodes2;
    numDimensions = size(nodes, 2);
    
    numKnots1_unique = size(unique(knots1), 1);
    numKnots2_unique = size(unique(knots2), 1);
    numKnots1_new = numKnots1 + q1 * numKnots1_unique;
    numKnots2_new = numKnots2 + q2 * numKnots2_unique;
    numNodes1_new = numNodes1 + q1 * (numKnots1_unique - 1);
    numNodes2_new = numNodes2 + q2 * (numKnots2_unique - 1);
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots1_new = knots1;
        knots2_new = knots2;
        nodes_new = nodes;
        p1_new = p1;
        p2_new = p2;
        
        return;
        
    elseif (q1 == 0 && q2 == 0)
        fprintf('Warning: Please specify a positive number for q1 or q2 for degree elevation.\n\nNo operation has been performed.\n\n');
        
        knots1_new = knots1;
        knots2_new = knots2;
        nodes_new = nodes;
        p1_new = p1;
        p2_new = p2;
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Degree elevation for direction 1
    %----------------------------------------------------------------------
    if (q1 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "line"
        %------------------------------------------------------------------
        temp2 = (0 : (numNodes2 - 1))';
        
        shifts_for_nodes     = numNodes1     * temp2;
        shifts_for_nodes_new = numNodes1_new * temp2;
        
        shifts_for_bezierNodes     = constant_p1p1     * temp2;
        shifts_for_bezierNodes_new = constant_p1pq1_p1 * temp2;
        
        clear temp2;
        
        
        %------------------------------------------------------------------
        %  Initialize the new knot vector and nodes array
        %------------------------------------------------------------------
        % Set xi1 to the first knot
        xi1 = knots1(1);
        
        knots1_new = zeros(numKnots1_new, 1);
        knots1_new((1 : constant_p1pq1_p1)', 1) = xi1;
        index_knot = constant_p1pq1_p1 + 1;
        
        nodes_new = zeros(numNodes1_new * numNodes2, numDimensions);
        nodes_new(1 + shifts_for_nodes_new, :) = nodes(1 + shifts_for_nodes, :);
        index_node = 2;
        
        
        %------------------------------------------------------------------
        %  Set variables for the Bezier element [xi1, xi2] x knots2
        %------------------------------------------------------------------
        % Set the Bezier nodes array to that for the first Bezier element
        bezierNodes = nodes(repmat((1 : constant_p1p1)', numNodes2, 1) + kron(shifts_for_nodes, ones(constant_p1p1, 1)), :);
        
        % Initialize the coefficients for knot insertion
        alphas_knotInsertion = zeros(p1 - 1, 1);
        
        % Initialize the Bezier nodes array for the degree elevated Bezier element
        bezierNodes_new = zeros(constant_p1pq1_p1 * numNodes2, numDimensions);
        
        % Bezier nodes shared by the next Bezier element
        bezierNodes_nextElement = zeros((p1 - 1) * numNodes2, numDimensions);
        
        % Indices of the knots corresponding to the first Bezier element
        a = constant_p1p1;
        b = a + 1;
        
        % Number of times we need to insert the knot xi2 to get a Bezier
        % element (this will be calculated later)
        numKnotsToInsert = 0;
        
        
        %------------------------------------------------------------------
        %  Compute the coefficients for the degree elevation
        %------------------------------------------------------------------
        % Array of coefficients used to degree elevate the Bezier elements
        alphas_degreeElevation = zeros(constant_p1pq1_p1, constant_p1p1);
        alphas_degreeElevation(1, 1) = 1;
        alphas_degreeElevation(constant_p1pq1_p1, constant_p1p1) = 1;
        
        % Compute the degree elevation coefficients
        for i = 2 : (floor(constant_p1pq1 / 2) + 1)
            temp = 1 / nchoosek(constant_p1pq1, i - 1);
            
            for j = max(1, i - q1) : min(constant_p1p1, i)
                alphas_degreeElevation(i, j) = temp * nchoosek(p1, j - 1) * nchoosek(q1, i - j);
            end
        end
        
        for i = (floor(constant_p1pq1 / 2) + 2) : constant_p1pq1
            for j = max(1, i - q1) : min(constant_p1p1, i)
                alphas_degreeElevation(i, j) = alphas_degreeElevation(constant_p1pq1_p1 - i + 1, constant_p1p1 - j + 1);
            end
        end
        
        
        %------------------------------------------------------------------
        %  Loop over the elements
        %------------------------------------------------------------------
        while (b < numKnots1)
            % Count how many times the knot xi2 is repeated
            i = b;
            
            while (b < numKnots1 && knots1(b) == knots1(b + 1))
                b = b + 1;
            end
            
            xi2 = knots1(b);
            
            numRepeats = b - i + 1;
            
            
            %--------------------------------------------------------------
            %  Knot insertion
            %--------------------------------------------------------------
            numKnotsToInsert_old = numKnotsToInsert;
            
            if (numKnotsToInsert_old > 0)
                lbz = 2 + floor(numKnotsToInsert_old / 2);
            else
                lbz = 2;
            end
            
            % Insert the knot xi2 to get a Bezier element
            numKnotsToInsert = p1 - numRepeats;
            
            if (numKnotsToInsert > 0)
                % Compute the coefficients for the knot insertion
                numerator = xi2 - xi1;
                
                for i = 1 : numKnotsToInsert
                    alphas_knotInsertion(i) = numerator / (knots1(a + numRepeats + i) - xi1);
                end
                
                % Compute the nodes for the Bezier element
                for i = 1 : numKnotsToInsert
                    s = numRepeats + i;
                    
                    for k = constant_p1p1 : -1 : (s + 1)
                        alpha = alphas_knotInsertion(k - s);
                        
                        bezierNodes(k + shifts_for_bezierNodes, :) = alpha * bezierNodes(k + shifts_for_bezierNodes, :) + (1 - alpha) * bezierNodes(k - 1 + shifts_for_bezierNodes, :);
                    end
                    
                    index = numKnotsToInsert - i + 1;
                    bezierNodes_nextElement(index + shifts_for_bezierNodes, :) = bezierNodes(constant_p1p1 + shifts_for_bezierNodes, :);
                end
                
                rbz = constant_p1pq1_p1 - floor((numKnotsToInsert + 1) / 2);
            else
                rbz = constant_p1pq1_p1;
            end
            
            
            %--------------------------------------------------------------
            %  Degree elevation
            %--------------------------------------------------------------
            % Compute the nodes for the degree elevated Bezier element
            for i = lbz : constant_p1pq1_p1
                bezierNodes_new(i + shifts_for_bezierNodes_new, :) = 0;
                
                for j = max(1, i - q1) : min(constant_p1p1, i)
                    bezierNodes_new(i + shifts_for_bezierNodes_new, :) = bezierNodes_new(i + shifts_for_bezierNodes_new, :) + alphas_degreeElevation(i, j) * bezierNodes(j + shifts_for_bezierNodes, :);
                end
            end
            
            
            %--------------------------------------------------------------
            %  Knot removal
            %--------------------------------------------------------------
            % Remove the knot xi1 (this was xi2 in the previous element)
            % from the new knot vector if there are more than needed to
            % meet the continuity requirement
            if (numKnotsToInsert_old > 1)
                first = index_knot - 2;
                last = index_knot;
                
                for t = 2 : numKnotsToInsert_old
                    i = first;
                    j = last;
                    k = j - index_knot + 2;
                    
                    while (j - i > t - 1)
                        % Update the nodes
                        if (i < index_node)
                            alpha = (xi2 - knots1_new(i)) / (xi1 - knots1_new(i));
                            
                            nodes_new(i + shifts_for_nodes_new, :) = alpha * nodes_new(i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(i - 1 + shifts_for_nodes_new, :);
                        end
                        
                        % Update the nodes for the next Bezier element
                        if (j >= lbz)
                            if (j - t + 1 <= index_knot - constant_p1pq1 + numKnotsToInsert_old)
                                alpha = (xi2 - knots1_new(j - t + 1)) / (xi2 - xi1);
                                
                                bezierNodes_new(k + shifts_for_bezierNodes_new, :) = alpha * bezierNodes_new(k + shifts_for_bezierNodes_new, :) + (1 - alpha) * bezierNodes_new(k + 1 + shifts_for_bezierNodes_new, :);
                            else
                                alpha = (xi2 - knots1_new(index_knot - 1)) / (xi2 - xi1);
                                
                                bezierNodes_new(k + shifts_for_bezierNodes_new, :) = alpha * bezierNodes_new(k + shifts_for_bezierNodes_new, :) + (1 - alpha) * bezierNodes_new(k + 1 + shifts_for_bezierNodes_new, :);
                            end
                        end
                        
                        i = i + 1;
                        j = j - 1;
                        k = k - 1;
                    end
                    
                    first = first - 1;
                    last = last + 1;
                end
            end
            
            
            %--------------------------------------------------------------
            %  Insert knot xi1 to the new knot vector
            %--------------------------------------------------------------
            if (a > constant_p1p1)
                knots1_new(index_knot + (0 : (constant_p1pq1 - numKnotsToInsert_old - 1))') = xi1;
                index_knot = index_knot + (constant_p1pq1 - numKnotsToInsert_old);
            end
            
            
            %--------------------------------------------------------------
            %  Insert the nodes to the new nodes array
            %--------------------------------------------------------------
            for i = lbz : rbz
                nodes_new(index_node + shifts_for_nodes_new, :) = bezierNodes_new(i + shifts_for_bezierNodes_new, :);
                index_node = index_node + 1;
            end
            
            
            %--------------------------------------------------------------
            %  Update indices for the next iteration
            %--------------------------------------------------------------
            if (b < numKnots1)
                for i = 1 : numKnotsToInsert
                    bezierNodes(i + shifts_for_bezierNodes, :) = bezierNodes_nextElement(i + shifts_for_bezierNodes, :);
                end
                
                for i = (numKnotsToInsert + 1) : constant_p1p1
                    bezierNodes(i + shifts_for_bezierNodes, :) = nodes(b - constant_p1p1 + i + shifts_for_nodes, :);
                end
                
                a = b;
                b = a + 1;
                
                xi1 = xi2;
            else
                % For the last element, insert knot xi2 to the new knot vector
                knots1_new(index_knot + (0 : constant_p1pq1)') = xi2;
                index_knot = index_knot + constant_p1pq1_p1;
            end
        end
        
        nodes = nodes_new;
        p1_new = constant_p1pq1;
        
        clear knots1 nodes_new;
        clear alphas_knotInsertion alphas_degreeElevation bezierNodes bezierNodes_new bezierNodes_nextElement;
        clear shifts_for_nodes shifts_for_nodes_new shifts_for_bezierNodes shifts_for_bezierNodes_new;
        
    % Default action for no degree elevation
    else
        knots1_new = knots1;
        p1_new = p1;
        
        clear knots1;
    end
    
    
    %----------------------------------------------------------------------
    %  Degree elevation for direction 2
    %----------------------------------------------------------------------
    if (q2 > 0)
        %------------------------------------------------------------------
        %  Shifts to find all the nodes on the same "line"
        %  
        %  Note that numNodes1 = numNodes1_new now. We add an extra shift
        %  of -(numNodes1_new - 1) so that we can sequentially traverse
        %  along direction 2 like we did along direction 1.
        %------------------------------------------------------------------
        shifts_for_nodes     = (0 : (numNodes1_new - 1))' - (numNodes1_new - 1);
        shifts_for_nodes_new = shifts_for_nodes;
        
        shifts_for_bezierNodes     = shifts_for_nodes;
        shifts_for_bezierNodes_new = shifts_for_nodes_new;
        
        
        %------------------------------------------------------------------
        %  Initialize the new knot vector and nodes array
        %------------------------------------------------------------------
        % Set xi1 to the first knot
        xi1 = knots2(1);
        
        knots2_new = zeros(numKnots2_new, 1);
        knots2_new((1 : constant_p2pq2_p1)', 1) = xi1;
        index_knot = constant_p2pq2_p1 + 1;
        
        nodes_new = zeros(numNodes1_new * numNodes2_new, numDimensions);
        nodes_new(numNodes1_new + shifts_for_nodes_new, :) = nodes(numNodes1_new + shifts_for_nodes, :);
        index_node = 2;
        
        
        %------------------------------------------------------------------
        %  Set variables for the Bezier element knots1 x [xi1, xi2]
        %------------------------------------------------------------------
        % Set the Bezier nodes array to that for the first Bezier element
        bezierNodes = nodes((1 : (numNodes1_new * constant_p2p1))', :);
        
        % Initialize the coefficients for knot insertion
        alphas_knotInsertion = zeros(p2 - 1, 1);
        
        % Initialize the Bezier nodes array for the degree elevated Bezier element
        bezierNodes_new = zeros(numNodes1_new * constant_p2pq2_p1, numDimensions);
        
        % Bezier nodes shared by the next Bezier element
        bezierNodes_nextElement = zeros(numNodes1_new * (p2 - 1), numDimensions);
        
        % Indices of the knots corresponding to the first Bezier element
        a = constant_p2p1;
        b = a + 1;
        
        % Number of times we need to insert the knot xi2 to get a Bezier
        % element (this will be calculated later)
        numKnotsToInsert = 0;
        
        
        %------------------------------------------------------------------
        %  Compute the coefficients for the degree elevation
        %------------------------------------------------------------------
        % Array of coefficients used to degree elevate the Bezier elements
        alphas_degreeElevation = zeros(constant_p2pq2_p1, constant_p2p1);
        alphas_degreeElevation(1, 1) = 1;
        alphas_degreeElevation(constant_p2pq2_p1, constant_p2p1) = 1;
        
        % Compute the degree elevation coefficients
        for i = 2 : (floor(constant_p2pq2 / 2) + 1)
            temp = 1 / nchoosek(constant_p2pq2, i - 1);
            
            for j = max(1, i - q2) : min(constant_p2p1, i)
                alphas_degreeElevation(i, j) = temp * nchoosek(p2, j - 1) * nchoosek(q2, i - j);
            end
        end
        
        for i = (floor(constant_p2pq2 / 2) + 2) : constant_p2pq2
            for j = max(1, i - q2) : min(constant_p2p1, i)
                alphas_degreeElevation(i, j) = alphas_degreeElevation(constant_p2pq2_p1 - i + 1, constant_p2p1 - j + 1);
            end
        end
        
        
        %------------------------------------------------------------------
        %  Loop over the elements
        %------------------------------------------------------------------
        while (b < numKnots2)
            % Count how many times the knot xi2 is repeated
            i = b;
            
            while (b < numKnots2 && knots2(b) == knots2(b + 1))
                b = b + 1;
            end
            
            xi2 = knots2(b);
            
            numRepeats = b - i + 1;
            
            
            %--------------------------------------------------------------
            %  Knot insertion
            %--------------------------------------------------------------
            numKnotsToInsert_old = numKnotsToInsert;
            
            if (numKnotsToInsert_old > 0)
                lbz = 2 + floor(numKnotsToInsert_old / 2);
            else
                lbz = 2;
            end
            
            % Insert the knot xi2 to get a Bezier element
            numKnotsToInsert = p2 - numRepeats;
            
            if (numKnotsToInsert > 0)
                % Compute the coefficients for the knot insertion
                numerator = xi2 - xi1;
                
                for i = 1 : numKnotsToInsert
                    alphas_knotInsertion(i) = numerator / (knots2(a + numRepeats + i) - xi1);
                end
                
                % Compute the nodes for the Bezier element
                for i = 1 : numKnotsToInsert
                    s = numRepeats + i;
                    
                    for k = constant_p2p1 : -1 : (s + 1)
                        alpha = alphas_knotInsertion(k - s);
                        
                        bezierNodes(numNodes1_new * k + shifts_for_bezierNodes, :) = alpha * bezierNodes(numNodes1_new * k + shifts_for_bezierNodes, :) + (1 - alpha) * bezierNodes(numNodes1_new * (k - 1) + shifts_for_bezierNodes, :);
                    end
                    
                    index = numKnotsToInsert - i + 1;
                    bezierNodes_nextElement(numNodes1_new * index + shifts_for_bezierNodes, :) = bezierNodes(numNodes1_new * constant_p2p1 + shifts_for_bezierNodes, :);
                end
                
                rbz = constant_p2pq2_p1 - floor((numKnotsToInsert + 1) / 2);
            else
                rbz = constant_p2pq2_p1;
            end
            
            
            %--------------------------------------------------------------
            %  Degree elevation
            %--------------------------------------------------------------
            % Compute the nodes for the degree elevated Bezier element
            for i = lbz : constant_p2pq2_p1
                bezierNodes_new(numNodes1_new * i + shifts_for_bezierNodes_new, :) = 0;
                
                for j = max(1, i - q2) : min(constant_p2p1, i)
                    bezierNodes_new(numNodes1_new * i + shifts_for_bezierNodes_new, :) = bezierNodes_new(numNodes1_new * i + shifts_for_bezierNodes_new, :) + alphas_degreeElevation(i, j) * bezierNodes(numNodes1_new * j + shifts_for_bezierNodes, :);
                end
            end
            
            
            %--------------------------------------------------------------
            %  Knot removal
            %--------------------------------------------------------------
            % Remove the knot xi1 (this was xi2 in the previous element)
            % from the new knot vector if there are more than needed to
            % meet the continuity requirement
            if (numKnotsToInsert_old > 1)
                first = index_knot - 2;
                last = index_knot;
                
                for t = 2 : numKnotsToInsert_old
                    i = first;
                    j = last;
                    k = j - index_knot + 2;
                    
                    while (j - i > t - 1)
                        % Update the nodes
                        if (i < index_node)
                            alpha = (xi2 - knots2_new(i)) / (xi1 - knots2_new(i));
                            
                            nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) = alpha * nodes_new(numNodes1_new * i + shifts_for_nodes_new, :) + (1 - alpha) * nodes_new(numNodes1_new * (i - 1) + shifts_for_nodes_new, :);
                        end
                        
                        % Update the nodes for the next Bezier element
                        if (j >= lbz)
                            if (j - t + 1 <= index_knot - constant_p2pq2 + numKnotsToInsert_old)
                                alpha = (xi2 - knots2_new(j - t + 1)) / (xi2 - xi1);
                                
                                bezierNodes_new(numNodes1_new * k + shifts_for_bezierNodes_new, :) = alpha * bezierNodes_new(numNodes1_new * k + shifts_for_bezierNodes_new, :) + (1 - alpha) * bezierNodes_new(numNodes1_new * (k - 1) + shifts_for_bezierNodes_new, :);
                            else
                                alpha = (xi2 - knots2_new(index_knot - 1)) / (xi2 - xi1);
                                
                                bezierNodes_new(numNodes1_new * k + shifts_for_bezierNodes_new, :) = alpha * bezierNodes_new(numNodes1_new * k + shifts_for_bezierNodes_new, :) + (1 - alpha) * bezierNodes_new(numNodes1_new * (k - 1) + shifts_for_bezierNodes_new, :);
                            end
                        end
                        
                        i = i + 1;
                        j = j - 1;
                        k = k - 1;
                    end
                    
                    first = first - 1;
                    last = last + 1;
                end
            end
            
            
            %--------------------------------------------------------------
            %  Insert knot xi1 to the new knot vector
            %--------------------------------------------------------------
            if (a > constant_p2p1)
                knots2_new(index_knot + (0 : (constant_p2pq2 - numKnotsToInsert_old - 1))') = xi1;
                index_knot = index_knot + (constant_p2pq2 - numKnotsToInsert_old);
            end
            
            
            %--------------------------------------------------------------
            %  Insert the nodes to the new nodes array
            %--------------------------------------------------------------
            for i = lbz : rbz
                nodes_new(numNodes1_new * index_node + shifts_for_nodes_new, :) = bezierNodes_new(numNodes1_new * i + shifts_for_bezierNodes_new, :);
                index_node = index_node + 1;
            end
            
            
            %--------------------------------------------------------------
            %  Update indices for the next iteration
            %--------------------------------------------------------------
            if (b < numKnots2)
                for i = 1 : numKnotsToInsert
                    bezierNodes(numNodes1_new * i + shifts_for_bezierNodes, :) = bezierNodes_nextElement(numNodes1_new * i + shifts_for_bezierNodes, :);
                end
                
                for i = (numKnotsToInsert + 1) : constant_p2p1
                    bezierNodes(numNodes1_new * i + shifts_for_bezierNodes, :) = nodes(numNodes1_new * (b - constant_p2p1 + i) + shifts_for_nodes, :);
                end
                
                a = b;
                b = a + 1;
                
                xi1 = xi2;
            else
                % For the last element, insert knot xi2 to the new knot vector
                knots2_new(index_knot + (0 : constant_p2pq2)') = xi2;
                index_knot = index_knot + constant_p2pq2_p1;
            end
        end
        
        p2_new = constant_p2pq2;
        
        clear knots2 nodes;
        clear alphas_knotInsertion alphas_degreeElevation bezierNodes bezierNodes_new bezierNodes_nextElement;
        clear shifts_for_nodes shifts_for_nodes_new shifts_for_bezierNodes shifts_for_bezierNodes_new;
        
    % Default action for no degree elevation
    else
        knots2_new = knots2;
        nodes_new = nodes;
        p2_new = p2;
        
        clear knots2 nodes;
    end
end