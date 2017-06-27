%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  This code is based on that by Mike Borden. (The code by Piegl and Tiller
%  in their NURBS book is incorrect and does not work for all cases.)
%  
%  
%  Summary:
%  
%  This routine does a p-refinement (degree elevation) of the given knot
%  vector in 1D.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open. The nodes are assumed to be
%  ordered along the direction of 1.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [knots_new, nodes_new, p_new] = refine_p_curve(knots, nodes, p, q);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      nodes is an array of nodes
%      p is the degree of the B-splines
%      q is the degree increment
%  
%  
%  Outputs:
%  
%  1. New knot vector
%  
%  2. New nodes array
%  
%  3. New degree
%--------------------------------------------------------------------------
function [knots_new, nodes_new, p_new] = refine_p_curve(knots, nodes, p, q)
    % Some useful constants
    constant_pp1 = p + 1;
    constant_ppq = p + q;
    constant_ppq_p1 = constant_ppq + 1;
    
    numKnots = size(knots, 1);
    numNodes = numKnots - constant_pp1;
    numDimensions = size(nodes, 2);
    
    numKnots_unique = size(unique(knots), 1);
    numKnots_new = numKnots + q * numKnots_unique;
    numNodes_new = numNodes + q * (numKnots_unique - 1);
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots_new = knots;
        nodes_new = nodes;
        p_new = p;
        
        return;
        
    elseif (q == 0)
        fprintf('Warning: Please specify a positive number for q for degree elevation.\n\nNo operation has been performed.\n\n');
        
        knots_new = knots;
        nodes_new = nodes;
        p_new = p;
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Degree elevation for direction 1
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %  Initialize the new knot vector and nodes array
    %----------------------------------------------------------------------
    % Set xi1 to the first knot
    xi1 = knots(1);
    
    knots_new = zeros(numKnots_new, 1);
    knots_new((1 : constant_ppq_p1)', 1) = xi1;
    index_knot = constant_ppq_p1 + 1;
    
    nodes_new = zeros(numNodes_new, numDimensions);
    nodes_new(1, :) = nodes(1, :);
    index_node = 2;
    
    
    %----------------------------------------------------------------------
    %  Set variables for the Bezier element [xi1, xi2]
    %----------------------------------------------------------------------
    % Set the Bezier nodes array to that for the first Bezier element
    bezierNodes = nodes((1 : constant_pp1)', :);
    
    % Initialize the coefficients for knot insertion
    alphas_knotInsertion = zeros(p - 1, 1);
    
    % Initialize the Bezier nodes array for the degree elevated Bezier element
    bezierNodes_new = zeros(constant_ppq_p1, numDimensions);
    
    % Bezier nodes shared by the next Bezier element
    bezierNodes_nextElement = zeros(p - 1, numDimensions);
    
    % Indices of the knots corresponding to the first Bezier element
    a = constant_pp1;
    b = a + 1;
    
    % Number of times we need to insert the knot xi2 to get a Bezier
    % element (this will be calculated later)
    numKnotsToInsert = 0;
    
    
    %----------------------------------------------------------------------
    %  Compute the coefficients for the degree elevation
    %----------------------------------------------------------------------
    % Array of coefficients used to degree elevate the Bezier elements
    alphas_degreeElevation = zeros(constant_ppq_p1, constant_pp1);
    alphas_degreeElevation(1, 1) = 1;
    alphas_degreeElevation(constant_ppq_p1, constant_pp1) = 1;
    
    % Compute the degree elevation coefficients
    for i = 2 : (floor(constant_ppq / 2) + 1)
        temp = 1 / nchoosek(constant_ppq, i - 1);
        
        for j = max(1, i - q) : min(constant_pp1, i)
            alphas_degreeElevation(i, j) = temp * nchoosek(p, j - 1) * nchoosek(q, i - j);
        end
    end
    
    for i = (floor(constant_ppq / 2) + 2) : constant_ppq
        for j = max(1, i - q) : min(constant_pp1, i)
            alphas_degreeElevation(i, j) = alphas_degreeElevation(constant_ppq_p1 - i + 1, constant_pp1 - j + 1);
        end
    end
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    while (b < numKnots)
        % Count how many times the knot xi2 is repeated
        i = b;
        
        while (b < numKnots && knots(b) == knots(b + 1))
            b = b + 1;
        end
        
        xi2 = knots(b);
        
        numRepeats = b - i + 1;
        
        
        %------------------------------------------------------------------
        %  Knot insertion
        %------------------------------------------------------------------
        numKnotsToInsert_old = numKnotsToInsert;
        
        if (numKnotsToInsert_old > 0)
            lbz = 2 + floor(numKnotsToInsert_old / 2);
        else
            lbz = 2;
        end
        
        % Insert the knot xi2 to get a Bezier element
        numKnotsToInsert = p - numRepeats;
        
        if (numKnotsToInsert > 0)
            % Compute the coefficients for the knot insertion
            numerator = xi2 - xi1;
            
            for i = 1 : numKnotsToInsert
                alphas_knotInsertion(i) = numerator / (knots(a + numRepeats + i) - xi1);
            end
            
            % Compute the nodes for the Bezier element
            for i = 1 : numKnotsToInsert
                s = numRepeats + i;
                
                for k = constant_pp1 : -1 : (s + 1)
                    alpha = alphas_knotInsertion(k - s);
                    
                    bezierNodes(k, :) = alpha * bezierNodes(k, :) + (1 - alpha) * bezierNodes(k - 1, :);
                end
                
                index = numKnotsToInsert - i + 1;
                bezierNodes_nextElement(index, :) = bezierNodes(constant_pp1, :);
            end
            
            rbz = constant_ppq_p1 - floor((numKnotsToInsert + 1) / 2);
        else
            rbz = constant_ppq_p1;
        end
        
        
        %------------------------------------------------------------------
        %  Degree elevation
        %------------------------------------------------------------------
        % Compute the nodes for the degree elevated Bezier element
        for i = lbz : constant_ppq_p1
            bezierNodes_new(i, :) = 0;
            
            for j = max(1, i - q) : min(constant_pp1, i)
                bezierNodes_new(i, :) = bezierNodes_new(i, :) + alphas_degreeElevation(i, j) * bezierNodes(j, :);
            end
        end
        
        
        %------------------------------------------------------------------
        %  Knot removal
        %------------------------------------------------------------------
        % Remove the knot xi1 (this was xi2 in the previous element) from
        % the new knot vector if there are more than needed to meet the
        % continuity requirement
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
                        alpha = (xi2 - knots_new(i)) / (xi1 - knots_new(i));
                        
                        nodes_new(i, :) = alpha * nodes_new(i, :) + (1 - alpha) * nodes_new(i - 1, :);
                    end
                    
                    % Update the nodes for the next Bezier element
                    if (j >= lbz)
                        if (j - t + 1 <= index_knot - constant_ppq + numKnotsToInsert_old)
                            alpha = (xi2 - knots_new(j - t + 1)) / (xi2 - xi1);
                            
                            bezierNodes_new(k, :) = alpha * bezierNodes_new(k, :) + (1 - alpha) * bezierNodes_new(k + 1, :);
                        else
                            alpha = (xi2 - knots_new(index_knot - 1)) / (xi2 - xi1);
                            
                            bezierNodes_new(k, :) = alpha * bezierNodes_new(k, :) + (1 - alpha) * bezierNodes_new(k + 1, :);
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
        
        
        %------------------------------------------------------------------
        %  Insert knot xi1 to the new knot vector
        %------------------------------------------------------------------
        if (a > constant_pp1)
            knots_new(index_knot + (0 : (constant_ppq - numKnotsToInsert_old - 1))') = xi1;
            index_knot = index_knot + (constant_ppq - numKnotsToInsert_old);
        end
        
        
        %------------------------------------------------------------------
        %  Insert the nodes to the new nodes array
        %------------------------------------------------------------------
        for i = lbz : rbz
            nodes_new(index_node, :) = bezierNodes_new(i, :);
            index_node = index_node + 1;
        end
        
        
        %------------------------------------------------------------------
        %  Update indices for the next iteration
        %------------------------------------------------------------------
        if (b < numKnots)
            for i = 1 : numKnotsToInsert
                bezierNodes(i, :) = bezierNodes_nextElement(i, :);
            end
            
            for i = (numKnotsToInsert + 1) : constant_pp1
                bezierNodes(i, :) = nodes(b - constant_pp1 + i, :);
            end
            
            a = b;
            b = a + 1;
            
            xi1 = xi2;
        else
            % For the last element, insert knot xi2 to the new knot vector
            knots_new(index_knot + (0 : constant_ppq)') = xi2;
            index_knot = index_knot + constant_ppq_p1;
        end
    end
    
    p_new = constant_ppq;
    
    clear knots nodes;
    clear alphas_knotInsertion alphas_degreeElevation bezierNodes bezierNodes_new bezierNodes_nextElement;
end