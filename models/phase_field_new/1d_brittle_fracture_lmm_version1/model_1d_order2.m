%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 1D cracked bar in tension, so that
%  we can validate the 2nd-order phase field theory for brittle fracture.
%  
%  We use Newton's method to solve for the displacement and phase fields.
%  In addition, Lagrangian multiplier method is used to try to enforce
%  the continuity requirement of the phase field at the crack.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_order2.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      restartTime is the Newton's iteration at which we start the simulation
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d_order2(path_to_assembly_directory, path_to_results_directory, restartTime)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  2nd-order phase field theory in 1D.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'         , ...
         'numNodesBeforePatch', ...
         'numMatrixEntries'   , ...
         'numDOFs'            , ...
         'numDOFsPerNode'     , ...
         ...
         'material_A'         , ...
         'material_G_c'       , ...
         'material_ell_0'     , ...
         ...
         'numConstraints'     , ...
         'maxNewtonsMethod'   , ...
         'tolNewtonsMethod');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'material_E'        , ...
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
         'LM_array1' , 'LM_array2' , 'LM_array3', ...
         'BCU_array1', 'BCU_array2', ...
         'BCF_array1', 'BCF_array2', ...
         'index_u'   , 'index_f');
    
    
    
    %----------------------------------------------------------------------
    %  Initialize the fields
    %----------------------------------------------------------------------
    if (restartTime == 0)
        initialize(2, path_to_assembly_directory, path_to_results_directory);
        load(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2', 'l0');
        
    else
        load(sprintf('%sfile_results_iteration%d', path_to_results_directory, restartTime), 'u1', 'u2', 'l0');
        
    end
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants for quadrature
    constant_p1p1 = p1 + 1;
    
    % Set the quadrature rule
    [z1, w] = set_1d_gauss_quadrature_for_bernstein(numQuadraturePoints);
    numQuadraturePointsPerElement = numQuadraturePoints;
    
    
    % Evaluate the Bernstein polynomials for direction 1 at the
    % quadrature points
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
    %   Begin: Loop over Newton's method (k)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_stiffness = material_E * material_A;
    constant_G_c       = material_G_c * material_A / (2 * material_ell_0);
    constant_ell_0_sq  = material_ell_0^2;
    
    
    % Total number of DOFs due to the displacement field, phase field,
    % and Lagrange multipliers
    numDOFs_total = (2 + numConstraints) * numDOFs;
    
    % Some useful constant for global assembly
    constant_p1p1_sq = constant_p1p1^2;
    
    
    % Initialization for element matrices
    elementMatrix = zeros(constant_p1p1);
    
    % Initialization for element vectors f1, f2
    elementVector = zeros(constant_p1p1, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(constant_p1p1, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : constant_p1p1_sq)';
    
    
    % Number of matrix entries that we compute for the tangent matrix
    numMatrixEntriesPerElement = (4 + 4 * numConstraints) * constant_p1p1_sq;
    
    % Initialize row, column, value arrays for the tangent matrix K
    rows_for_K    = zeros(numMatrixEntriesPerElement * numElements1, 1);
    columns_for_K = zeros(numMatrixEntriesPerElement * numElements1, 1);
    values_for_K  = zeros(numMatrixEntriesPerElement * numElements1, 1);
    
    
    % Flag that indicates whether we are close to convergence
    flag_isCloseToConvergence = 0;
    
    % Norm of the residual vector
    residual_norm = 0;
    
    
    for k = restartTime : maxNewtonsMethod
        % Initialize the internal force
        f = zeros(numDOFs_total, 1);
        
        % External force due to increments in body force and tractions
        f(BCF_array1(:, 1)          ) = -BCF_array1(:, 2);
        f(BCF_array2(:, 1) + numDOFs) = -BCF_array2(:, 2);
        
        % Index that was last used to set the entry of K
        lastIndex_for_K = 0;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
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
            %   Begin: Loop over quadrature points (q, q_e)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Initialize the element tangent matrices and the element
            % internal force vectors
            K11_e = elementMatrix;
            K12_e = elementMatrix;
            K13_e = elementMatrix;
            
            K22_e = elementMatrix;
            K23_e = elementMatrix;
            
            f1_e  = elementVector;
            f2_e  = elementVector;
            f3_e  = elementVector;
            
            
            % Find the equation indices for the element
            index_equation1 = LM_array1(:, e);
            index_equation2 = LM_array2(:, e);
            index_equation3 = LM_array3(:, e);
            
            % Get the current Newton's guess for the displacement field,
            % phase field, and Lagrange multipliers
            u1_e = u1(index_equation1)';
            u2_e = u2(index_equation2)';
            l0_e = l0(index_equation3)';
            
            
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
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u_der1 = u1_e * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                c_der0 = u2_e * basis_physical_der0;
                c_der1 = u2_e * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the Lagrange multiplier field
                %----------------------------------------------------------
                l0_der0 = l0_e * basis_physical_der0;
                
                
                %----------------------------------------------------------
                %  Form the element matrix K11
                %----------------------------------------------------------
                % The subscripts denote the derivatives of N_A and N_B
                C11 = constant_stiffness * c_der0^2;
                
                % Additional terms from Lagrange multiplier method
                C11 = C11 + (2 * l0_der0 * c_der0);
                
                K11_e = K11_e + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K12
                %----------------------------------------------------------
                % The subscripts denote the derivatives of N_A and N_B
                C10 = 2 * constant_stiffness * c_der0 * u_der1;
                
                % Additional terms from Lagrange multiplier method
                C10 = C10 + (2 * l0_der0 * u_der1);
                
                K12_e = K12_e + (w(q_e) * Jacobian * C10) * (basis_physical_der1 * basis_physical_der0');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K13
                %----------------------------------------------------------
                % Additional terms from Lagrange multiplier method
                C10 = 2 * c_der0 * u_der1;
                
                K13_e = K13_e + (w(q_e) * Jacobian * C10) * (basis_physical_der1 * basis_physical_der0');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K22
                %----------------------------------------------------------
                % The subscripts denote the derivatives of N_A and N_B
                C00 = constant_stiffness * u_der1^2 + constant_G_c;
                C11 = constant_G_c * (4 * constant_ell_0_sq);
                
                K22_e = K22_e + (w(q_e) * Jacobian * C00) * (basis_physical_der0 * basis_physical_der0') ...
                              + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K23
                %----------------------------------------------------------
                % Additional terms from Lagrange multiplier method
                C00 = u_der1^2;
                
                K23_e = K23_e + (w(q_e) * Jacobian * C00) * (basis_physical_der0 * basis_physical_der0');
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f1
                %----------------------------------------------------------
                % The subscript denotes the derivative of N_A
                C1 = constant_stiffness * c_der0^2 * u_der1;
                
                % Additional terms from Lagrange multiplier method
                C1 = C1 + (2 * l0_der0 * c_der0 * u_der1);
                
                f1_e = f1_e + (w(q_e) * Jacobian * C1) * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f2
                %----------------------------------------------------------
                % The subscript denotes the derivative of N_A
                C0 = constant_stiffness * c_der0 * u_der1^2 - constant_G_c * (1 - c_der0);
                C1 = constant_G_c * (4 * constant_ell_0_sq * c_der1);
                
                % Additional terms from Lagrange multiplier method
                C0 = C0 + (l0_der0 * u_der1^2);
                
                f2_e = f2_e + (w(q_e) * Jacobian * C0) * basis_physical_der0 ...
                            + (w(q_e) * Jacobian * C1) * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f3
                %----------------------------------------------------------
                % Additional terms from Lagrange multiplier method
                C0 = c_der0 * u_der1^2;
                
                f3_e = f3_e + (w(q_e) * Jacobian * C0) * basis_physical_der0;
                
            end
            
            
            % Do we set the off-diagonal blocks to zero?
            K12_e = elementMatrix;
%            K13_e = elementMatrix;
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Global assembly
            %--------------------------------------------------------------
            % Shift the indices
            index_equation2 = index_equation2 +     numDOFs;
            index_equation3 = index_equation3 + 2 * numDOFs;
            
            
            % Add the element matrix K11
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation1);
            columns_for_K(index) = kron(index_equation1, vector_of_ones);
            values_for_K(index)  = reshape(K11_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K12
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation1);
            columns_for_K(index) = kron(index_equation2, vector_of_ones);
            values_for_K(index)  = reshape(K12_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K13
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation1);
            columns_for_K(index) = kron(index_equation3, vector_of_ones);
            values_for_K(index)  = reshape(K13_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K21
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation2);
            columns_for_K(index) = kron(index_equation1, vector_of_ones);
            values_for_K(index)  = reshape(K12_e', constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K22
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation2);
            columns_for_K(index) = kron(index_equation2, vector_of_ones);
            values_for_K(index)  = reshape(K22_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K23
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation2);
            columns_for_K(index) = kron(index_equation3, vector_of_ones);
            values_for_K(index)  = reshape(K23_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K31
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation3);
            columns_for_K(index) = kron(index_equation1, vector_of_ones);
            values_for_K(index)  = reshape(K13_e', constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K32
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K(index)    = kron(vector_of_ones, index_equation3);
            columns_for_K(index) = kron(index_equation2, vector_of_ones);
            values_for_K(index)  = reshape(K23_e', constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element matrix K33
            % Do nothing (K33 = 0)
            
            
            % Add the element RHS vector f1
            f(index_equation1) = f(index_equation1) + f1_e;
            
            
            % Add the element RHS vector f2
            f(index_equation2) = f(index_equation2) + f2_e;
            
            
            % Add the element RHS vector f3
            f(index_equation3) = f(index_equation3) + f3_e;
            
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Check for convergence
        %  
        %  We accept that the Newton's guess corresponds to an equilibrium
        %  state via a two-step check.
        %------------------------------------------------------------------
        if (k > 0)
            % Save the coefficients for debugging
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, k), 'u1', 'u2', 'l0', '-v6');
            
            
            fprintf('\n');
            fprintf('  Close to convergence = %d\n', flag_isCloseToConvergence);
            fprintf('  Newton''s increment abs., 2-norm = %.4e\n', increment_norm);
            fprintf('  Residual vector    abs., 2-norm = %.4e\n', residual_norm);
            
            
            % Check whether the current Newton's increment is small,
            % and if so, whether the current residual is smaller than
            % the previous residual
            if (flag_isCloseToConvergence == 1 && residual_norm < residual_norm0)
                fprintf('\n');
                fprintf('  Newton''s method converged at iteration %d.\n', k);
                
                break;
                
            elseif (k == maxNewtonsMethod)
                fprintf('\n');
                fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                
                break;
                
            end
        end
        
        
        %------------------------------------------------------------------
        %  Solve for the Newton's increment
        %------------------------------------------------------------------
        % Assemble the tangent matrix
        K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs_total, numDOFs_total);
        
        
        % Save the norm of the residual vector at the previous iteration
        residual_norm0 = residual_norm;
        
        % Evaluate the norm of the residual vector at the current iteration
        residual_norm = norm(f(index_u), 2);
        
        
        % Apply a row preconditioner
%        precond = spdiags(1./max(K, [], 2), 0, numDOFs_total, numDOFs_total);
%        K = precond * K;
%        f = precond * f;
        
%        clear precond;
        
        
        % Solve for the Newton's increment
        u_increment = zeros(numDOFs_total, 1);
        u_increment(index_u) = -K(index_u, index_u) \ f(index_u);
        
        clear K f;
        
        
        % Update the Newton's guess
        u1 = u1 + u_increment((1 : numDOFs)'              );
        u2 = u2 + u_increment((1 : numDOFs)' +     numDOFs);
        l0 = l0 + u_increment((1 : numDOFs)' + 2 * numDOFs);
        
        
        % Evaluate the norm of the Newton's increment and check if
        % we are close to convergence
        increment_norm = norm(u_increment(index_u), 2);
        
        if (increment_norm < tolNewtonsMethod)
            flag_isCloseToConvergence = 1;
        end
        
        clear u_increment;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over Newton's method
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
