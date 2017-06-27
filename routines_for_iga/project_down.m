%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine finds the nodes so that the B-spline in (d + 1)-dimensional
%  space can be projected to the NURBS in the d-dimensional space.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      nodes = project_down(nodes_projectUp);
%  
%  where,
%  
%      nodes_projectUp is an array of B-spline nodes
%  
%  
%  Output:
%  
%  1. Array of NURBS nodes
%--------------------------------------------------------------------------
function nodes = project_down(nodes_projectUp)
    numNodes = size(nodes_projectUp, 1);
    numDimensions = size(nodes_projectUp, 2) - 1;
    
    % Initialize the array of B-spline nodes
    nodes = zeros(numNodes, numDimensions + 1);
    
    for j = 1 : numDimensions
        % Compute the new coordinates
        nodes(:, j) = nodes_projectUp(:, j) ./ nodes_projectUp(:, numDimensions + 1);
    end
    
    % Save the weights
    nodes(:, numDimensions + 1) = nodes_projectUp(:, numDimensions + 1);
end