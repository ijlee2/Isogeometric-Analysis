%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the phase field equation for 1D cracked bar in
%  tension.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_order8.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d_order8(path_to_assembly_directory, path_to_results_directory)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  8th-order phase field theory in 1D.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'         , ...
         'numNodesBeforePatch', ...
         'numMatrixEntries'   , ...
         'numDOFs'            , ...
         'numDOFsPerNode'     , ...
         ...
         'material_ell_0');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'elementType'       , ...
         'p1'                , ...
         'nodes'             , ...
         'numNodesPerElement', ...
         'numDOFsPerElement' , ...
         'bezierExtractions1', ...
         'numElements1'      , ...
         'elementSizes1'     , ...
         'IEN_array'         , ...
         'numQuadraturePoints');
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'LM_array2'      , ...
         'BCU_array2'     , ...
         'BCF_array2'     , ...
         'numUnknownDOFs2', ...
         'index_u2'       , ...
         'index_f2');
    
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants for quadrature
    constant_p1p1 = p1 + 1;
    
    % Set the quadrature rule
    [z1, w] = set_1d_gauss_quadrature_for_bernstein(numQuadraturePoints);
    numQuadraturePointsPerElement = numQuadraturePoints;
    
    
    % Evaluate the Bernstein polynomials for direction 1 at the
    % quadrature points, and at the two endpoints
    z1 = [z1; 0; 1];
    temp = zeros(constant_p1p1, numQuadraturePoints + 2);
    
    Bernstein1_der0 = temp;
    Bernstein1_der1 = temp;
    Bernstein1_der2 = temp;
    Bernstein1_der3 = temp;
    Bernstein1_der4 = temp;
    
    for q_e = 1 : numQuadraturePoints + 2
        temp = eval_1d_bernstein_der(z1(q_e), p1);
        
        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
        Bernstein1_der2(:, q_e) = temp(:, 3);
        Bernstein1_der3(:, q_e) = temp(:, 4);
        Bernstein1_der4(:, q_e) = temp(:, 5);
    end
    
    clear z1 temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over alterations
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_ell_0_sq = material_ell_0^2;
    
    
    % Some useful constant for global assembly
    constant_p1p1_sq = constant_p1p1^2;
    
    
    % Modification for Lagrange multiplier
    numConstraints = p1;
    numDOFs_total = numDOFs + numConstraints;
    index_u2 = [index_u2; (numDOFs + 1 : numDOFs + numConstraints)'];
    
    % Initialization for solution and RHS vectors
    globalVector  = zeros(numDOFs_total, 1);
    
    % Initialization for element matrices K11, K12, K22
    elementMatrix = zeros(constant_p1p1);
    
    % Initialization for element vectors f1, f2
    elementVector = zeros(constant_p1p1, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(constant_p1p1, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : constant_p1p1_sq)';
    
    
    % Number of matrix entries that we compute for the momentum and
    % phase field equations
    numMatrixEntriesPerElement = constant_p1p1_sq;
    
    % Initialize row, column, value arrays for a stiffness matrix K
    temp = zeros(numMatrixEntriesPerElement * numElements1 + 2 * numConstraints * constant_p1p1, 1);
    
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
    
    clear temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Find the phase field (u2)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Initialize the RHS vector and prescribe the known forces
    f2 = globalVector;
    f2(BCF_array2(:, 1)) = BCF_array2(:, 2);
    
    % Index that was last used to set the entry of K
    lastIndex_for_K = 0;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over elements (e = e1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for e = 1 : numElements1
        % Find the Bezier extraction matrix
        bezierExtractions1_e = bezierExtractions1(:, :, e);
        
        % Evaluate the map derivative dt/dxi (constant)
        dt_dxi = 1 / elementSizes1(e);
        
        
        %------------------------------------------------------------------
        %  Evaluate the basis functions in the parametric domain
        %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
        %------------------------------------------------------------------
        % Find the positions of the nodes
        nodes_e = nodes(IEN_array(:, e), :)';
        
        % Evaluate the B-splines
        basis_der0 =             bezierExtractions1_e  * Bernstein1_der0;
        basis_der1 = (dt_dxi   * bezierExtractions1_e) * Bernstein1_der1;
        basis_der2 = (dt_dxi^2 * bezierExtractions1_e) * Bernstein1_der2;
        basis_der3 = (dt_dxi^3 * bezierExtractions1_e) * Bernstein1_der3;
        basis_der4 = (dt_dxi^4 * bezierExtractions1_e) * Bernstein1_der4;
        
        
        %------------------------------------------------------------------
        %  Evaluate the map derivatives (x = x1)
        %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
        %------------------------------------------------------------------
        dx_dxi = nodes_e * basis_der1;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over quadrature points (q_e)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Initialize the element stiffness matrix and the element
        % internal force vector
        K22_e = elementMatrix;
        f2_e  = elementVector;
        
        % Find the equation indices for the element
        index_equation2 = LM_array2(:, e);
        
        
        for q_e = 1 : numQuadraturePointsPerElement
            %--------------------------------------------------------------
            %  Form the Jacobian matrix
            %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
            %--------------------------------------------------------------
            JacobianMatrix = dx_dxi(:, q_e);
            
            % Evaluate the Jacobian
            Jacobian = JacobianMatrix / dt_dxi;
            
            if (Jacobian <= 0)
                fprintf('\n');
                fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                
                quit;
            end
            
            
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the physical domain
            %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
            %--------------------------------------------------------------
            JacobianMatrix_inv = 1 / JacobianMatrix;
            
            basis_physical_der0 =                        basis_der0(:, q_e);
            basis_physical_der1 = JacobianMatrix_inv   * basis_der1(:, q_e);
            basis_physical_der2 = JacobianMatrix_inv^2 * basis_der2(:, q_e);
            basis_physical_der3 = JacobianMatrix_inv^3 * basis_der3(:, q_e);
            basis_physical_der4 = JacobianMatrix_inv^4 * basis_der4(:, q_e);
            
            
            %--------------------------------------------------------------
            %  Form the element matrix K22
            %--------------------------------------------------------------
            % Coefficients for the (regular) phase field theory.
            % The subscripts denote the derivatives of N_A and N_B.
            C00 = 1;
            C11 =         constant_ell_0_sq;
            C22 = 3/8   * constant_ell_0_sq^2;
            C33 = 1/16  * constant_ell_0_sq^3;
            C44 = 1/256 * constant_ell_0_sq^4;
            
            K22_e = K22_e + (w(q_e) * Jacobian * C00) * (basis_physical_der0 * basis_physical_der0') ...
                          + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1') ...
                          + (w(q_e) * Jacobian * C22) * (basis_physical_der2 * basis_physical_der2') ...
                          + (w(q_e) * Jacobian * C33) * (basis_physical_der3 * basis_physical_der3') ...
                          + (w(q_e) * Jacobian * C44) * (basis_physical_der4 * basis_physical_der4');
            
            
            %--------------------------------------------------------------
            %  Form the element RHS vector f2
            %--------------------------------------------------------------
            % Coefficients for the (regular) phase field theory.
            % The subscript denotes the derivative of N_A.
            C0 = 1;
            
            f2_e = f2_e + (w(q_e) * Jacobian * C0) * basis_physical_der0;
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over quadrature points
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Add the element matrix K22
        index = lastIndex_for_K + indices_for_elementMatrix;
        
        rows_for_K(index)    = kron(vector_of_ones, index_equation2);
        columns_for_K(index) = kron(index_equation2, vector_of_ones);
        values_for_K(index)  = reshape(K22_e, constant_p1p1_sq, 1);
        
        lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
        
        
        % Add the element RHS vector f2
        f2(index_equation2) = f2(index_equation2) + f2_e;
        
        
        % Evaluate the B-splines at the endpoints
        x_left  = nodes_e * basis_der0(:, numQuadraturePoints + 1);
        x_right = nodes_e * basis_der0(:, numQuadraturePoints + 2);
        
        % Check if the crack is located inside the element
        if (x_left <= 0 && 0 < x_right)
            if (abs(x_left) < 1e-12)
                % Form the Jacobian matrix
                JacobianMatrix = dx_dxi(:, numQuadraturePoints + 1);
                
                for m = 0 : 3
                    % Evaluate the basis functions in the physical domain
                    JacobianMatrix_inv = 1 / JacobianMatrix;
                    
                    if (m == 0)
                        basis_physical_derm =                        basis_der0(:, numQuadraturePoints + 1);
                        
                    elseif (m == 1)
                        basis_physical_derm = JacobianMatrix_inv   * basis_der1(:, numQuadraturePoints + 1);
                        
                    elseif (m == 2)
                        basis_physical_derm = JacobianMatrix_inv^2 * basis_der2(:, numQuadraturePoints + 1);
                        
                    elseif (m == 3)
                        basis_physical_derm = JacobianMatrix_inv^3 * basis_der3(:, numQuadraturePoints + 1);
                        
                    end
                    
                    
                    % Add the element matrices K23 due to Lagrange multiplier
                    index = lastIndex_for_K + (1 : constant_p1p1)';
                    
                    rows_for_K(index)    = index_equation2;
                    columns_for_K(index) = (numDOFs + m + 1) * ones(constant_p1p1, 1);
                    values_for_K(index)  = basis_physical_derm;
                    
                    lastIndex_for_K = lastIndex_for_K + constant_p1p1;
                    
                    
                    % Add the element matrices K32 due to Lagrange multiplier
                    index = lastIndex_for_K + (1 : constant_p1p1)';
                    
                    rows_for_K(index)    = (numDOFs + m + 1) * ones(constant_p1p1, 1);
                    columns_for_K(index) = index_equation2;
                    values_for_K(index)  = basis_physical_derm;
                    
                    lastIndex_for_K = lastIndex_for_K + constant_p1p1;
                    
                end
                
            elseif (abs(x_right) < 1e-12)
                % Form the Jacobian matrix
                JacobianMatrix = dx_dxi(:, numQuadraturePoints + 2);
                
                for m = 0 : 3
                    % Evaluate the basis functions in the physical domain
                    JacobianMatrix_inv = 1 / JacobianMatrix;
                    
                    if (m == 0)
                        basis_physical_derm =                        basis_der0(:, numQuadraturePoints + 2);
                        
                    elseif (m == 1)
                        basis_physical_derm = JacobianMatrix_inv   * basis_der1(:, numQuadraturePoints + 2);
                        
                    elseif (m == 2)
                        basis_physical_derm = JacobianMatrix_inv^2 * basis_der2(:, numQuadraturePoints + 2);
                        
                    elseif (m == 3)
                        basis_physical_derm = JacobianMatrix_inv^3 * basis_der3(:, numQuadraturePoints + 2);
                        
                    end
                    
                    
                    % Add the element matrices K23 due to Lagrange multiplier
                    index = lastIndex_for_K + (1 : constant_p1p1)';
                    
                    rows_for_K(index)    = index_equation2;
                    columns_for_K(index) = (numDOFs + m + 1) * ones(constant_p1p1, 1);
                    values_for_K(index)  = basis_physical_derm;
                    
                    lastIndex_for_K = lastIndex_for_K + constant_p1p1;
                    
                    
                    % Add the element matrices K32 due to Lagrange multiplier
                    index = lastIndex_for_K + (1 : constant_p1p1)';
                    
                    rows_for_K(index)    = (numDOFs + m + 1) * ones(constant_p1p1, 1);
                    columns_for_K(index) = index_equation2;
                    values_for_K(index)  = basis_physical_derm;
                    
                    lastIndex_for_K = lastIndex_for_K + constant_p1p1;
                    
                end
                
            else
                fprintf('\nError: x = 0 is not located on either endpoint.\n\n');
                
                quit;
                
            end
        end
        
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Solve for the phase field
    %----------------------------------------------------------------------
    % Assemble the stiffness matrix
    K22 = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs_total, numDOFs_total);
    
    
    % Apply a row preconditioner
    precond = spdiags(1./max(K22, [], 2), 0, numDOFs_total, numDOFs_total);
    K22 = precond * K22;
    f2  = precond * f2;
    
    
    % Prescribe the known coefficients and solve for the unknown
    u2 = globalVector;
    u2(BCU_array2(:, 1)) = BCU_array2(:, 2);
    u2(index_u2) = K22(index_u2, index_u2) \ (f2(index_u2) - K22(index_u2, index_f2) * u2(index_f2));
    
    clear K22 f2;
    
    
    fprintf('\n');
    fprintf('  Equation 2 (phase field equation) has been solved.\n\n');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Save the results
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    save(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u2', '-v6');
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
