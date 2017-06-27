%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine draws the NURBS curve that is specified by the knot vector
%  of degree p and the nodes (control points) in 1D, 2D, or 3D.
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
%      draw_nurbs_curve(knots, nodes, p);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      nodes is an array of nodes
%      p is the degree of the NURBS
%  
%  
%  Output:
%  
%  1. Plot of the NURBS curve
%--------------------------------------------------------------------------
function draw_nurbs_curve(knots, nodes, p)
    % Some useful constants
    constant_pp1 = p + 1;
    
    numKnots = size(knots, 1);
    numNodes = numKnots - constant_pp1;
    numDimensions = size(nodes, 2) - 1;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vector, there should be %d nodes.\n\n', numNodes);
        
        return;
        
    elseif (numDimensions < 1 || numDimensions > 3)
        fprintf('Error: The nodes must be in the 1D, 2D, or 3D space.\n\n');
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extraction matrices
    %----------------------------------------------------------------------
    [bezierExtractions, nodeIndexShifts, numElements] = build_bezier_extraction(knots, p);
    
    
    %----------------------------------------------------------------------
    %  Build the IEN array
    %----------------------------------------------------------------------
    IEN_array = build_ien_array(knots, [], [], p, [], [], nodeIndexShifts, [], []);
    
    
    %----------------------------------------------------------------------
    %  Project NURBS to the B-spline in the higher dimension
    %----------------------------------------------------------------------
    nodes_projectUp = project_up(nodes);
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    figure;
    
    % Number of points where the fields are evaluated
    if (numElements > 8)
        numPointsPerElementCoordinate = 9;
    else
        numPointsPerElementCoordinate = 17;
    end
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElementCoordinate)';
    
    % Evaluate the Bernstein polynomials
    Bernstein_der0 = zeros(constant_pp1, numPointsPerElementCoordinate);
    
    for j = 1 : numPointsPerElementCoordinate
        Bernstein_der0(:, j) = eval_1d_bernstein(t(j), p);
    end
    
    % Initialize the points on the physical domain
    x = zeros(numPointsPerElementCoordinate, 1);
    y = zeros(numPointsPerElementCoordinate, 1);
    if (numDimensions == 3)
        z = zeros(numPointsPerElementCoordinate, 1);
    end
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Find the positions of the nodes
        nodes_e = nodes_projectUp(IEN_array(:, e), :)';
        
        % Find the Bezier extraction matrix
        bezierExtraction_e = bezierExtractions(:, :, e);
        
        % Evaluate the B-splines
        Bspline_der0 = bezierExtraction_e * Bernstein_der0;
        
        
        %------------------------------------------------------------------
        %  Loop over the points
        %------------------------------------------------------------------
        for j = 1 : numPointsPerElementCoordinate
            % Evaluate the basis functions for the physical domain
            basis_der0 = Bspline_der0(:, j);
            
            % Find the weight
            w = nodes_e(numDimensions + 1, :) * basis_der0;
            
            % Find the corresponding point on the physical domain
            x(j) = (nodes_e(1, :) * basis_der0) / w;
            
            if (numDimensions >= 2)
                y(j) = (nodes_e(2, :) * basis_der0) / w;
            end
            
            if (numDimensions == 3)
                z(j) = (nodes_e(3, :) * basis_der0) / w;
            end
        end
        
        
        %------------------------------------------------------------------
        %  Draw the NURBS curve
        %------------------------------------------------------------------
        if (numDimensions <= 2)
            %{
            % Uncomment this block to draw the control mesh
            if (numDimensions == 2)
                nodes_e = nodes(IEN_array(:, e), :);
                
                plot(nodes_e(:, 1), nodes_e(:, 2), '--o', ...
                     'LineWidth', 2, 'Color', [0.75; 0.75; 0.75], 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.80; 0.80; 0.80]); hold on;
            end
            %}
            
            % Draw the element
            plot(x, y, '-', 'LineWidth', 2, 'Color', [0.70; 0.35; 0.45]); hold on;
            
            % Draw the end points
            plot(x(1), y(1), 's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
            plot(x(numPointsPerElementCoordinate), y(numPointsPerElementCoordinate), 's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
            
        else
            %{
            % Uncomment this block to draw the control mesh
            nodes_e = nodes(IEN_array(:, e), :);
            
            plot3(nodes_e(:, 1), nodes_e(:, 2), nodes_e(:, 3), '--o', ...
                  'LineWidth', 2, 'Color', [0.75; 0.75; 0.75], 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.80; 0.80; 0.80]); hold on;
            %}
            
            % Draw the element
            plot3(x, y, z, '-', 'LineWidth', 2, 'Color', [0.70; 0.35; 0.45]); hold on;
            
            % Draw the end points
            plot3(x(1), y(1), z(1), 's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
            plot3(x(numPointsPerElementCoordinate), y(numPointsPerElementCoordinate), z(numPointsPerElementCoordinate), 's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
            
        end
    end
    
    xlabel('x', 'FontSize', 90);
    ylabel('y', 'FontSize', 90);
    zlabel('z', 'FontSize', 90, 'Rotation', 0);
    axis square;
    grid on;
    set(gca, 'FontSize', 54);
    if (numDimensions <= 2)
        ylabel('y', 'FontSize', 90, 'Rotation', 0);
        view(2);
    end
end