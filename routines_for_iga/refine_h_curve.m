%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine does a h-refinement (knot refinement) of the given knot
%  vector in 1D.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open, and the knots for insertion are
%  assumed to be in the interior of the knot vector. The nodes are assumed
%  to be ordered along the direction of 1.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [knots_new, nodes_new] = refine_h_curve(knots, nodes, p, knotsForInsertion);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      nodes is an array of nodes
%      p is the degree of the B-splines
%      knotsForInsertion is the knot vector that we seek to insert
%  
%  
%  Outputs:
%  
%  1. New knot vector (same degree)
%  
%  2. New nodes array
%--------------------------------------------------------------------------
function [knots_new, nodes_new] = refine_h_curve(knots, nodes, p, knotsForInsertion)
    % Some useful constants
    constant_pp1 = p + 1;
    
    numKnots = size(knots, 1);
    numNodes = numKnots - constant_pp1;
    
    numKnotInsertions = size(knotsForInsertion, 1);
    numKnots_new = numKnots + numKnotInsertions;
%   numNodes_new = numNodes + numKnotInsertions;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        knots_new = knots;
        nodes_new = nodes;
        
        return;
        
    elseif (numKnotInsertions == 0)
        fprintf('Warning: Please specify at least one knot for insertion.\n\nNo operation has been performed.\n\n');
        
        knots_new = knots;
        nodes_new = nodes;
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Knot refinement for direction 1
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %  Initialize the new knot vector and nodes array
    %----------------------------------------------------------------------
    knots_new = [knots; zeros(numKnotInsertions, 1)];
    
    clear knots;
    
    nodes_new = [nodes; repmat(nodes(numNodes, :), numKnotInsertions, 1)];
    
    clear nodes;
    
    
    %----------------------------------------------------------------------
    %  Loop over the knots for insertion
    %----------------------------------------------------------------------
    % Search where the knot xi belongs at this index
    knotSearch_begin = constant_pp1;
    
    for b = 1 : numKnotInsertions
        % Knot that we will insert
        xi = knotsForInsertion(b);
        
        % Find where to insert the knot
        for k = knotSearch_begin : (numKnots_new - 1)
            if (knots_new(k) <= xi && xi < knots_new(k + 1))
                index_knot = k + 1;
                
                break;
            end
        end
        
        % Shift the "end" nodes
        for i = (numNodes + b - 1) : -1 : index_knot
            nodes_new(i, :) = nodes_new(i - 1, :);
        end
        
        % Update the "middle" nodes
        for i = (index_knot - 1) : -1 : (index_knot - p)
            alpha = (xi - knots_new(i)) / (knots_new(i + p) - knots_new(i));
            
            nodes_new(i, :) = alpha * nodes_new(i, :) + (1 - alpha) * nodes_new(i - 1, :);
        end
        
        % Shift the "end" knots
        for i = (numKnots + b) : -1 : (index_knot + 1)
            knots_new(i) = knots_new(i - 1);
        end
        
        
        %------------------------------------------------------------------
        %  Insert knot xi to the new knot vector
        %------------------------------------------------------------------
        knots_new(index_knot) = xi;
        
        
        %------------------------------------------------------------------
        %  Update index for the next iteration
        %------------------------------------------------------------------
        knotSearch_begin = index_knot;
    end
end