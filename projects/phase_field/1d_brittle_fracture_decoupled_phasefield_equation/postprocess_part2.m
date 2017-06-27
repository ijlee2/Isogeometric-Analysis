%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the surface energy for the decoupled 1D bar
%  problem.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess_part2.m output
%  
%  Note that,
%  
%      order is the order of the phase field theory (2, 4, 6, 8)
%      path_to_assembly_directory is the path to the assembly directory
%      path_to_results_directory is the path to the results directory
%      path_to_outputs_directory is the path to the plots directory
%  
%  
%  Output:
%  
%  1. Errors in the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function postprocess_part2(order, path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    %----------------------------------------------------------------------
    %  Initialize the arrays
    %----------------------------------------------------------------------
    % Set the refinement levels
    numRefinements = [0; 1; 2; 3; 4];
    
    temp = zeros(size(numRefinements, 1), 1);
    
    surface_energy = temp;
    
    clear temp;
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants for quadrature
    p1 = order / 2;
    constant_p1p1 = p1 + 1;
    numQuadraturePoints = constant_p1p1;
    
    % Set the quadrature rule
    [z1, w] = set_1d_gauss_quadrature_for_bernstein(numQuadraturePoints);
    numQuadraturePointsPerElement = numQuadraturePoints;
    
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der2 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der3 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der4 = zeros(constant_p1p1, numQuadraturePoints);
    
    for q_e = 1 : numQuadraturePoints
        temp = eval_1d_bernstein_der(z1(q_e), p1);

        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
        if (order >= 4)
            Bernstein1_der2(:, q_e) = temp(:, 3);
        end
        if (order >= 6)
            Bernstein1_der3(:, q_e) = temp(:, 4);
        end
        if (order >= 8)
            Bernstein1_der4(:, q_e) = temp(:, 5);
        end
    end
    
    clear z1 temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Compute the errors and energies
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for i = 1 : size(numRefinements, 1)
        % Load the global assembly file
        load(sprintf('%snumRefinements%d/file_assembly_global', path_to_assembly_directory, numRefinements(i)), ...
             'material_ell_0');
        
        % Load the patch assembly file
        load(sprintf('%snumRefinements%d/file_assembly_patch%d', path_to_assembly_directory, numRefinements(i), 1), ...
             'nodes'             , ...
             'bezierExtractions1', ...
             'numElements1'      , ...
             'elementSizes1'     , ...
             'IEN_array');
        
        % Load the BCs
        load(sprintf('%snumRefinements%d/file_bc', path_to_assembly_directory, numRefinements(i)), 'LM_array2');
        
        % Load the results file
        load(sprintf('%snumRefinements%d/file_results_alternation%d', path_to_results_directory, numRefinements(i), 0), 'u2');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Some useful constants for phase field theory
        constant_G_c      = 1 / (4 * material_ell_0);
        constant_ell_0_sq = material_ell_0^2;
        
        
        for e = 1 : numElements1
            % Find the Bezier extraction matrix
            bezierExtractions1_e = bezierExtractions1(:, :, e);
            
            % Evaluate the map derivative dt/dxi (constant)
            dt_dxi = 1 / elementSizes1(e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the parametric domain
            %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            % Find the positions of the nodes
            nodes_e = nodes(IEN_array(:, e), :)';
            
            % Evaluate the B-splines
            basis_der0     =             bezierExtractions1_e  * Bernstein1_der0;
            basis_der1     = (dt_dxi   * bezierExtractions1_e) * Bernstein1_der1;
            if (order >= 4)
                basis_der2 = (dt_dxi^2 * bezierExtractions1_e) * Bernstein1_der2;
            end
            if (order >= 6)
                basis_der3 = (dt_dxi^3 * bezierExtractions1_e) * Bernstein1_der3;
            end
            if (order >= 8)
                basis_der4 = (dt_dxi^4 * bezierExtractions1_e) * Bernstein1_der4;
            end
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx_dxi = nodes_e * basis_der1;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over quadrature points (q_e, q)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Find the coefficients of the solution vectors
            u2_e = u2(LM_array2(:, e))';
            
            % Initialize the surface energy in the element
            surface_energy_e = 0;
            
            
            for q_e = 1 : numQuadraturePointsPerElement
                %----------------------------------------------------------
                %  Form the Jacobian matrix
                %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                %----------------------------------------------------------
                JacobianMatrix = dx_dxi(:, q_e);
                
                % Evaluate the Jacobian
                Jacobian = JacobianMatrix / dt_dxi;
                
                if (Jacobian <= 0)
                    fprintf('\n');
                    fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                    
                    quit;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the physical domain
                %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                %----------------------------------------------------------
                JacobianMatrix_inv = 1 / JacobianMatrix;
                
                basis_physical_der0     =                        basis_der0(:, q_e);
                basis_physical_der1     = JacobianMatrix_inv   * basis_der1(:, q_e);
                if (order >= 4)
                    basis_physical_der2 = JacobianMatrix_inv^2 * basis_der2(:, q_e);
                end
                if (order >= 6)
                    basis_physical_der3 = JacobianMatrix_inv^3 * basis_der3(:, q_e);
                end
                if (order >= 8)
                    basis_physical_der4 = JacobianMatrix_inv^4 * basis_der4(:, q_e);
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                u2_der0     = u2_e * basis_physical_der0;
                u2_der1     = u2_e * basis_physical_der1;
                if (order >= 4)
                    u2_der2 = u2_e * basis_physical_der2;
                end
                if (order >= 6)
                    u2_der3 = u2_e * basis_physical_der3;
                end
                if (order >= 8)
                    u2_der4 = u2_e * basis_physical_der4;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the surface energy
                %----------------------------------------------------------
                switch order
                    case 2
                        crack_density = (1 - u2_der0)^2 + 4   * constant_ell_0_sq * u2_der1^2;
                        
                    case 4
                        crack_density = (1 - u2_der0)^2 + 2   * constant_ell_0_sq * u2_der1^2 +         constant_ell_0_sq^2 * u2_der2^2;
                        
                    case 6
                        crack_density = (1 - u2_der0)^2 + 4/3 * constant_ell_0_sq * u2_der1^2 + 16/27 * constant_ell_0_sq^2 * u2_der2^2 + 64/729 * constant_ell_0_sq^3 * u2_der3^2;
                        
                    case 8
                        crack_density = (1 - u2_der0)^2 +       constant_ell_0_sq * u2_der1^2 + 3/8   * constant_ell_0_sq^2 * u2_der2^2 + 1/16   * constant_ell_0_sq^3 * u2_der3^2 + 1/256 * constant_ell_0_sq^4 * u2_der4^2;
                        
                end
                
                surface_energy_e = surface_energy_e + (w(q_e) * Jacobian * (constant_G_c * crack_density));
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Add the element contribution
            %--------------------------------------------------------------
            surface_energy(i) = surface_energy(i) + surface_energy_e;
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Compute the errors and energies
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Save the results
    %----------------------------------------------------------------------
    save(sprintf('%ssurface_energy', path_to_outputs_directory), 'surface_energy', '-v7.3');
end
