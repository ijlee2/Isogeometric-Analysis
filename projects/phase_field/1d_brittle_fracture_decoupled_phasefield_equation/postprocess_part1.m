%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots the phase field for the decoupled 1D bar problem.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess_part1.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      path_to_outputs_directory is the path to the plots directory
%  
%  
%  Output:
%  
%  1. Plots of the displacement and phase fields (.png files)
%--------------------------------------------------------------------------
function postprocess_part1(path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory, numRefinements)
    close all;
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'material_ell_0');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'p1'                , ...
         'nodes'             , ...
         'bezierExtractions1', ...
         'numElements1'      , ...
         'elementSizes1'     , ...
         'IEN_array'         , ...
         'numQuadraturePoints');
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), 'LM_array2');
    
    % Load the results file
    load(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0));
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants
    constant_p1p1 = p1 + 1;
    
    
    % Number of points that we consider on each element
    numPointsPerElement = 11;
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElement)';
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numPointsPerElement);

    for q_e = 1 : numPointsPerElement
        Bernstein1_der0(:, q_e) = eval_1d_bernstein(t(q_e), p1);
    end

    clear t;
    
    
    % Specify the elements to zoom in on
    elementIndex_begin = numElements1 / 2 - (200 * 2^numRefinements);
    elementIndex_end   = numElements1 / 2 + (200 * 2^numRefinements);
    numElementsToConsider = elementIndex_end - elementIndex_begin + 1;
    
    % Initialize the points on the physical domain
    temp = zeros(numPointsPerElement * numElementsToConsider, 1);
    
    x       = temp;
    u2_der0 = temp;
    
    clear temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over elements (e = e1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Counter for the point
    q = 1;
    
    
    for e = 1 : numElements1
        % Find the Bezier extraction matrix
        bezierExtractions1_e = bezierExtractions1(:, :, e);
        
        
        %------------------------------------------------------------------
        %  Evaluate the basis functions in the parametric domain
        %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
        %------------------------------------------------------------------
        % Find the positions of the nodes
        nodes_e = nodes(IEN_array(:, e), :)';
        
        % Evaluate the B-splines
        basis_der0 = bezierExtractions1_e * Bernstein1_der0;
        
        
        %------------------------------------------------------------------
        %  Loop over points
        %------------------------------------------------------------------
        % Find the coefficients of the solution vectors
        u2_e = u2(LM_array2(:, e))';
        
        
        for q_e = 1 : numPointsPerElement
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the physical domain
            %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
            %--------------------------------------------------------------
            basis_physical_der0 = basis_der0(:, q_e);
            
            
            %--------------------------------------------------------------
            %  Find the point on the physical domain
            %--------------------------------------------------------------
            x(q) = nodes_e * basis_physical_der0;
            
            
            %--------------------------------------------------------------
            %  Evaluate the displacement and phase fields
            %--------------------------------------------------------------
            u2_der0(q) = u2_e * basis_physical_der0;
            
            
            % Increment the counter
            q = q + 1;
        end
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Plot the exact solution
    %----------------------------------------------------------------------
    u2_der0_exact = function_c_exact(x, material_ell_0, 2 * p1);
    h(2) = plot(x, u2_der0_exact, '--', 'Color', [0.00, 0.00, 0.00], 'LineWidth', 3); hold on;
    
    
    %----------------------------------------------------------------------
    %  Plot the phase field
    %----------------------------------------------------------------------
    h(1) = plot(x, u2_der0, '-', 'Color', [0.70, 0.35, 0.45], 'LineWidth', 3);
    xlabel('x', 'FontSize', 90);
    ylabel('c(x)', 'FontSize', 90);
    xmin = -1;
    xmax = 1;
    axis([xmin xmax -0.05 1.05]);
    grid on;
    
    legend(h, {' FE solution', ' analytical solution'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, 1, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf('%splot_c_alternation%d.png', path_to_outputs_directory, 0));
end


function c = function_c_exact(x, ell_0, order)
    % Normalize the physical domain
    y = abs(x)/ell_0;
    
    switch order
        case 2
            c = 1 - exp(-1/2*y);
            
        case 4
            c = 1 - exp(-y) .* (1 + y);
            
        case 6
            c = 1 - exp(-3/2*y) .* (1 + 3/2*y + 9/8*y.^2);
            
        case 8
            c = 1 - exp(-2*y) .* (1 + 2*y + 2*y.^2 + 4/3*y.^3);
            
    end
end
