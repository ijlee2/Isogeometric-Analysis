%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots the displacement and phase fields at each alternation
%  step for the problem of 1D cracked bar in tension.
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
function postprocess_part1(path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'material_L'    , ...
         'material_A'    , ...
         'material_G_c'  , ...
         'material_ell_0', ...
         ...
         'maxAlternations');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'material_E'        , ...
         'p1'                , ...
         'nodes'             , ...
         'bezierExtractions1', ...
         'numElements1'      , ...
         'elementSizes1'     , ...
         'IEN_array'         , ...
         'numQuadraturePoints');
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), 'LM_array1', 'LM_array2', 'u_L');
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants
    constant_p1p1 = p1 + 1;
    
    
    % Number of points that we consider on each element
    numPointsPerElement = 5;
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElement)';
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numPointsPerElement);

    for q_e = 1 : numPointsPerElement
        Bernstein1_der0(:, q_e) = eval_1d_bernstein(t(q_e), p1);
    end

    clear t;
    
    
    % Initialize the points on the physical domain
    temp = zeros(numPointsPerElement * numElements1, 1);
    
    x       = temp;
    u1_der0 = temp;
    u2_der0 = temp;
    
    clear temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Plot the fields
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Counter for figures
    count_figure = 1;
    
    
    for alternation = 0 : maxAlternations
        % Load the results file
        load(sprintf('%sfile_results_alternation%d', path_to_results_directory, alternation));
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counter for the point
        q = 1;
        
        
        for e = 1 : numElements1
            % Find the Bezier extraction matrix
            bezierExtractions1_e = bezierExtractions1(:, :, e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the parametric domain
            %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            % Find the positions of the nodes
            nodes_e = nodes(IEN_array(:, e), :)';
            
            % Evaluate the B-splines
            basis_der0 = bezierExtractions1_e * Bernstein1_der0;
            
            
            %--------------------------------------------------------------
            %  Loop over points
            %--------------------------------------------------------------
            % Find the coefficients of the solution vectors
            u1_e = u1(LM_array1(:, e))';
            u2_e = u2(LM_array2(:, e))';
            
            
            for q_e = 1 : numPointsPerElement
                %----------------------------------------------------------
                %  Evaluate the basis functions in the physical domain
                %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                %----------------------------------------------------------
                basis_physical_der0 = basis_der0(:, q_e);
                
                
                %----------------------------------------------------------
                %  Find the point on the physical domain
                %----------------------------------------------------------
                if (alternation == 0)
                    x(q) = nodes_e * basis_physical_der0;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement and phase fields
                %----------------------------------------------------------
                u1_der0(q) = u1_e * basis_physical_der0;
                u2_der0(q) = u2_e * basis_physical_der0;
                
                
                % Increment the counter
                q = q + 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Plot the displacement field
        %------------------------------------------------------------------
        figure(count_figure);
        
        plot(x, u1_der0, '-', 'LineWidth', 3);
        xlabel('x', 'FontSize', 90);
        ylabel('u(x)', 'FontSize', 90);
        xmin = -material_L;
        xmax = material_L;
        axis([xmin xmax -0.05*u_L 1.05*u_L]);
        grid on;
        
        set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, u_L, 11)');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
        print('-dpng', '-r300', sprintf('%splot_u_alternation%d.png', path_to_outputs_directory, alternation));
        
        % Increment the counter for figure
        count_figure = count_figure + 1;
        
        
        %------------------------------------------------------------------
        %  Plot the phase field
        %------------------------------------------------------------------
        figure(count_figure);
        
        plot(x, u2_der0, '-', 'LineWidth', 3);
        xlabel('x', 'FontSize', 90);
        ylabel('c(x)', 'FontSize', 90);
        xmin = -material_L;
        xmax = material_L;
        axis([xmin xmax -0.05 1.05]);
        grid on;
        
        set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, 1, 11)');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
        print('-dpng', '-r300', sprintf('%splot_c_alternation%d.png', path_to_outputs_directory, alternation));
        
        % Increment the counter for figure
        count_figure = count_figure + 1;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Plot the fields
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
end
