%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes (1) the L2 errors between the analytical and FE
%  solutions of the displacement and phase fields, and (2) the strain,
%  surface, and potential energies at the last alternation step.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess_part3.m output
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
function postprocess_part3(order, p, path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    %----------------------------------------------------------------------
    %  Initialize the arrays
    %----------------------------------------------------------------------
    % Set the refinement levels
    numRefinements = [0; 1; 2; 3; 4];
    
    temp = zeros(size(numRefinements, 1), 1);
    
    strain_energy = temp;
    
    clear temp;
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants for quadrature
    p1 = p;
    constant_p1p1 = p1 + 1;
    numQuadraturePoints = constant_p1p1;
    
    % Set the quadrature rule
    [z1, w] = set_1d_gauss_quadrature_for_bernstein(numQuadraturePoints);
    numQuadraturePointsPerElement = numQuadraturePoints;
    
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints);
    
    for q_e = 1 : numQuadraturePoints
        temp = eval_1d_bernstein_der(z1(q_e), p1);

        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
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
             'material_L', ...
             'material_A', ...
             'material_G_c', ...
             'material_ell_0', ...
             ...
             'maxAlternations');
        
        % Load the patch assembly file
        load(sprintf('%snumRefinements%d/file_assembly_patch%d', path_to_assembly_directory, numRefinements(i), 1), ...
             'material_E'        , ...
             'nodes'             , ...
             'bezierExtractions1', ...
             'numElements1'      , ...
             'elementSizes1'     , ...
             'IEN_array');
        
        % Load the BCs
        load(sprintf('%snumRefinements%d/file_bc', path_to_assembly_directory, numRefinements(i)), 'LM_array1', 'LM_array2', 'u_L');
        
        % Load the results file
        load(sprintf('%snumRefinements%d/file_results_alternation%d', path_to_results_directory, numRefinements(i), 0));
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Some useful constants
        constant_stiffness = material_E * material_A;
        
        
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
            basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
            basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
            
            
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
            u1_e = u1(LM_array1(:, e))';
            
            % Initialize the errors and energies in the element
            strain_energy_e = 0;
            
            
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
                
                basis_physical_der0 =                      basis_der0(:, q_e);
                basis_physical_der1 = JacobianMatrix_inv * basis_der1(:, q_e);
                
                
                %----------------------------------------------------------
                %  Find the point on the physical domain
                %----------------------------------------------------------
                x = nodes_e * basis_physical_der0;
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u1_der1 = u1_e * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                y = abs(x) / material_ell_0;
                
                switch order
                    case 2
                        c_der0 = 1 - exp(-1/2*y);
                        
                    case 4
                        c_der0 = 1 - exp(-y) .* (1 + y);
                        
                    case 6
                        c_der0 = 1 - exp(-3/2*y) .* (1 + 3/2*y + 9/8*y.^2);
                        
                    case 8
                        c_der0 = 1 - exp(-2*y) .* (1 + 2*y + 2*y.^2 + 4/3*y.^3);
                        
                end

                
                %----------------------------------------------------------
                %  Evaluate the strain energy
                %----------------------------------------------------------
                strain_energy_e = strain_energy_e + (w(q_e) * Jacobian * (c_der0^2 * 0.5 * constant_stiffness * u1_der1^2));
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Add the element contribution
            %--------------------------------------------------------------
            strain_energy(i) = strain_energy(i) + strain_energy_e;
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
    save(sprintf('%serrors_and_energies', path_to_outputs_directory), 'numRefinements', 'strain_energy', '-v7.3');
end
