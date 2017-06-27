%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine finds the nodes so that the NURBS in d-dimensional space
%  can be projected to the B-spline in the (d + 1)-dimensional space.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      nodes_projectUp = project_up(nodes);
%  
%  where,
%  
%      nodes is an array of NURBS nodes
%  
%  
%  Output:
%  
%  1. Array of B-spline nodes
%--------------------------------------------------------------------------
function nodes_projectUp = project_up(nodes)
    numNodes = size(nodes, 1);
    numDimensions = size(nodes, 2) - 1;
    
    % Initialize the array of B-spline nodes
    nodes_projectUp = zeros(numNodes, numDimensions + 1);
    
    for j = 1 : numDimensions
        % Compute the new coordinates
        nodes_projectUp(:, j) = nodes(:, j) .* nodes(:, numDimensions + 1);
    end
    
    % Save the weights
    nodes_projectUp(:, numDimensions + 1) = nodes(:, numDimensions + 1);
end