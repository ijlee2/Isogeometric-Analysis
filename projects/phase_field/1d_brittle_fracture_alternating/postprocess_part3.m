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
function postprocess_part3(order, path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    %----------------------------------------------------------------------
    %  Initialize the arrays
    %----------------------------------------------------------------------
    % Set the refinement levels
    numRefinements = [0; 1; 2; 3; 4];
    
    temp = zeros(size(numRefinements, 1), 1);
    
    L2_error_u       = temp;
    L2_error_c       = temp;
    strain_energy    = temp;
    surface_energy   = temp;
    potential_energy = temp;
    
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
        load(sprintf('%snumRefinements%d/file_results_alternation%d', path_to_results_directory, numRefinements(i), maxAlternations), 'u1', 'u2');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Some useful constants
        constant_stiffness = material_E * material_A;
        constant_G_c       = material_G_c * material_A / (4 * material_ell_0);
        constant_ell_0_sq  = material_ell_0^2;
        
        
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
            u1_e = u1(LM_array1(:, e))';
            u2_e = u2(LM_array2(:, e))';
            
            % Initialize the errors and energies in the element
            L2_error_u_e       = 0;
            L2_error_c_e       = 0;
            strain_energy_e    = 0;
            surface_energy_e   = 0;
            potential_energy_e = 0;
            
            
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
                %  Find the point on the physical domain
                %----------------------------------------------------------
                x = nodes_e * basis_physical_der0;
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u1_der0 = u1_e * basis_physical_der0;
                u1_der1 = u1_e * basis_physical_der1;
                
                
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
                %  Evaluate the analytical solutions
                %----------------------------------------------------------
                u_exact = function_u_exact(x, u_L);
                c_exact = function_c_exact(x, material_ell_0, order);
                
                
                
                %----------------------------------------------------------
                %  Evaluate the L2 error of the displacement field
                %----------------------------------------------------------
                L2_error_u_e = L2_error_u_e + (w(q_e) * Jacobian * (u_exact - u1_der0)^2);
                
                
                %----------------------------------------------------------
                %  Evaluate the L2 error of the phase field
                %----------------------------------------------------------
                L2_error_c_e = L2_error_c_e + (w(q_e) * Jacobian * (c_exact - u2_der0)^2);
                
                
                %----------------------------------------------------------
                %  Evaluate the strain energy
                %----------------------------------------------------------
                degradation_function = u2_der0^2;
                
                strain_energy_e = strain_energy_e + (w(q_e) * Jacobian * (degradation_function * 0.5 * constant_stiffness * u1_der1^2));
                
                
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
                
                
                %----------------------------------------------------------
                %  Evaluate the potential energy
                %----------------------------------------------------------
                potential_energy_e = strain_energy_e + surface_energy_e;
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Add the element contribution
            %--------------------------------------------------------------
            L2_error_u(i)       = L2_error_u(i)       + L2_error_u_e      ;
            L2_error_c(i)       = L2_error_c(i)       + L2_error_c_e      ;
            strain_energy(i)    = strain_energy(i)    + strain_energy_e   ;
            surface_energy(i)   = surface_energy(i)   + surface_energy_e  ;
            potential_energy(i) = potential_energy(i) + potential_energy_e;
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
    end
    
    
    % Take the square root for the L2 errors
    L2_error_u = sqrt(L2_error_u);
    L2_error_c = sqrt(L2_error_c);
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Compute the errors and energies
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Save the results
    %----------------------------------------------------------------------
    save(sprintf('%serrors_and_energies', path_to_outputs_directory), 'numRefinements', 'L2_error_u', 'L2_error_c', 'strain_energy', 'surface_energy', 'potential_energy', '-v7.3');
end


function u = function_u_exact(x, u_L)
    if (x <= 0)
        u = 0;
    else
        u = u_L;
    end
end


function c = function_c_exact(x, ell_0, order)
    % Normalize the physical domain
    y = abs(x)/ell_0;
    
    switch order
        case 2
            c = 1 - exp(-1/2*y);
            
        case 4
            c = 1 - exp(-y) * (1 + y);
            
        case 6
            c = 1 - exp(-3/2*y) * (1 + 3/2*y + 9/8*y^2);
            
        case 8
            c = 1 - exp(-2*y) * (1 + 2*y + 2*y^2 + 4/3*y^3);
            
    end
end
