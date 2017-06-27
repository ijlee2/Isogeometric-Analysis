%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 1D cracked bar in tension, so that
%  we can validate the 2nd-order phase field theory for brittle fracture.
%  The momentum and phase field equations are solved alternatingly.
%  
%  In addition, augmented Lagrangian method is used to try to enforce the
%  continuity requirement at the crack.
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
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d_order2(path_to_assembly_directory, path_to_results_directory)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  2nd-order phase field theory in 1D.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'              , ...
         'numNodesBeforePatch'     , ...
         'numMatrixEntries'        , ...
         'numDOFs'                 , ...
         'numDOFsPerNode'          , ...
         'GN_array'                , ...
         ...
         'material_A'              , ...
         'material_G_c'            , ...
         'material_ell_0'          , ...
         'maxAlternations'         , ...
         ...
         'maxNewtonsMethod'        , ...
         'tolNewtonsMethod_u'      , ...
         'tolNewtonsMethod_c'      , ...
         'LMM_beta');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'elementType'       , ...
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
         'LM_array1', 'LM_array2'  , ...
         'BCU_array1', 'BCU_array2', ...
         'BCF_array1', 'BCF_array2', ...
         'index_u1', 'index_u2'    , ...
         'index_f1', 'index_f2');
    
    
    
    %----------------------------------------------------------------------
    %  Initialize the fields
    %----------------------------------------------------------------------
    initialize(2, path_to_assembly_directory, path_to_results_directory);
    
    load(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda0');
    
    
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
    %   Begin: Loop over alterations
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_stiffness = material_E * material_A;
    constant_G_c = material_G_c * material_A / (2 * material_ell_0);
    constant_ell_0_sq = material_ell_0^2;
    
    
    % Row, column, value arrays for the tangent matrices K11 and K22
    temp = zeros(numMatrixEntries, 1);
    
    rows_for_K11    = temp;
    columns_for_K11 = temp;
    values_for_K11  = temp;
    
    rows_for_K22    = temp;
    columns_for_K22 = temp;
    values_for_K22  = temp;
    
    clear temp;
    
    
    % Number of matrix entries that we compute per element
    numMatrixEntriesPerElement = numDOFsPerElement^2;
    
    % Vector of zeros for Newton's method
    vector_of_zeros = zeros(numDOFs, 1);
    
    % Vector of ones for global assembly
    vector_of_ones = ones(numDOFsPerElement, 1);
    
    
    for alternation = 1 : maxAlternations
        fprintf('\n');
        fprintf('- Alternation index = %d\n', alternation);
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   1. Find the displacement field (u1)
        %   
        %   Begin: Loop over Newton's method (k)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % By definition, the initial Newton's guess is the solution vector
        % from the previous alternation step
        u1_new = u1;
        
        % Initialize the Newton's increment
        u1_increment = vector_of_zeros;
        
        % External force vector due to body force and tractions
        f1_external = vector_of_zeros;
        f1_external(BCF_array1(:, 1)) = BCF_array1(:, 2);
        
        
        % Flag that indicates whether we are close to convergence
        flag_isCloseToConvergence = 0;
        
        % Norm of the residual vector
        residual_norm = 0;
        
        
        for k = 0 : maxNewtonsMethod
            % Initialize the internal force
            f1_internal = vector_of_zeros;
            
            % Initialize the indices of the matrix entries
            index_for_K11 = (1 : numMatrixEntriesPerElement)';
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over elements (e = e1)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Counters for the quadrature point
            q = 1;
            
            
            for e = 1 : numElements1
                % Find the Bezier extraction matrix
                bezierExtractions1_e = bezierExtractions1(:, :, e);
                
                % Evaluate the map derivative dt/dxi (constant)
                dt_dxi = 1 / elementSizes1(e);
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the parametric domain
                %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                % Find the positions of the nodes
                nodes_e = nodes(IEN_array(:, e), :)';
                
                % Evaluate the B-splines
                basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
                basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the map derivatives (x = x1)
                %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                dx_dxi = nodes_e * basis_der1;
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over quadrature points (q, q_e)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Initialize the element tangent matrix and the element
                % internal force vector
                K11_e = zeros(numDOFsPerElement);
                f1_e  = zeros(numDOFsPerElement, 1);
                
                % Find the equation indices for the element
                index_equation1 = LM_array1(:, e);
                index_equation2 = LM_array2(:, e);
                
                % Get the current Newton's guess for the solution vector
                u1_e = u1_new(index_equation1)';
                u2_e = u2(index_equation2)';
                
                
                for q_e = 1 : numQuadraturePointsPerElement
                    %------------------------------------------------------
                    %  Form the Jacobian matrix
                    %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                    %------------------------------------------------------
                    JacobianMatrix = dx_dxi(:, q_e);
                    
                    % Evaluate the Jacobian
                    Jacobian = JacobianMatrix / dt_dxi;
                    
                    if (Jacobian <= 0)
                        fprintf('\n');
                        fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                        
                        quit;
                    end
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the basis functions in the physical domain
                    %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                    %------------------------------------------------------
                    JacobianMatrix_inv = 1 / JacobianMatrix;
                    
                    basis_physical_der0 =                      basis_der0(:, q_e);
                    basis_physical_der1 = JacobianMatrix_inv * basis_der1(:, q_e);
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the displacement field
                    %------------------------------------------------------
                    u_der1 = u1_e * basis_physical_der1;
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the phase field
                    %------------------------------------------------------
                    c_der0 = u2_e * basis_physical_der0;
                    
                    
                    %------------------------------------------------------
                    %  Form the element matrix K11
                    %------------------------------------------------------
                    % Some useful constants for augmented Lagrangian method
                    constant0 = constant_stiffness * c_der0^2;
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscripts denote the derivatives of N_A and N_B.
                    C11 = constant0;
                    
                    
                    % Additional terms from the augmented Lagrangian method
                    C11 = C11 + (lambda0(q) * constant0) + LMM_beta * (3/2 * constant0^2 * u_der1^2);
                    
                    
                    K11_e = K11_e + w(q_e) * Jacobian * (C11 * (basis_physical_der1 * basis_physical_der1'));
                    
                    
                    %------------------------------------------------------
                    %  Form the element RHS vector f1
                    %------------------------------------------------------
                    % Some useful constants for augmented Lagrangian method
                    constant0 = constant_stiffness * c_der0^2 * u_der1;
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscript denotes the derivative of N_A.
                    C1 = constant0;
                    
                    
                    % Additional terms from the augmented Lagrangian method
                    C1 = C1 + (lambda0(q) * constant0) + LMM_beta * (1/2 * constant0^2 * u_der1);
                    
                    
                    f1_e = f1_e + w(q_e) * Jacobian * (C1 * basis_physical_der1);
                    
                    
                    % Increment the counter for quadrature point
                    q = q + 1;
                end
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Loop over quadrature points
                % ---------------------------------------------------------
                %----------------------------------------------------------
                
                
                %----------------------------------------------------------
                %  Global assembly
                %----------------------------------------------------------
                % Add the element matrix K11
                rows_for_K11(index_for_K11)    = kron(vector_of_ones, index_equation1);
                columns_for_K11(index_for_K11) = kron(index_equation1, vector_of_ones);
                values_for_K11(index_for_K11)  = reshape(K11_e, numMatrixEntriesPerElement, 1);
                index_for_K11 = index_for_K11 + numMatrixEntriesPerElement;
                
                
                % Add the element RHS vector f1
                f1_internal(index_equation1) = f1_internal(index_equation1) + f1_e;
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over elements
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Check for convergence
            %  
            %  We accept that the Newton's guess corresponds to an
            %  equilibrium state via a two-step check.
            %--------------------------------------------------------------
            if (k > 0)
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
                    fprintf('  Equation 1 (momentum equation) has been solved.\n\n');
                    
                    % Update the displacement field
                    u1 = u1_new;
                    
                    break;
                    
                elseif (k == maxNewtonsMethod)
                    fprintf('\n');
                    fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                    fprintf('  Equation 1 (momentum equation) has been solved.\n\n');
                    
                    % Update the displacement field
                    u1 = u1_new;
                    
                    break;
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Solve for the Newton's increment
            %--------------------------------------------------------------
            % Assemble the tangent matrix
            K11 = sparse(rows_for_K11, columns_for_K11, values_for_K11, numDOFs, numDOFs);
            
            
            % Save the norm of the residual vector at the previous iteration
            residual_norm0 = residual_norm;
            
            % Form the residual vector (note the negative sign)
            residual = f1_external - f1_internal;
            
            % Evaluate the norm of the residual vector at the current iteration
            residual_norm = norm(residual(index_u1), 2);
            
            
            % Apply a row preconditioner
            precond = spdiags(1./max(K11, [], 2), 0, numDOFs, numDOFs);
            K11 = precond * K11;
            residual = precond * residual;
            
            
            % Upon initialization, we told the solution vectors to satisfy
            % the displacements prescribed on the boundary. Hence, we tell
            % the Newton's increment to satisfy zero displacements on the
            % boundary.
            u1_increment(index_u1) = K11(index_u1, index_u1) \ residual(index_u1);
            
            
            % Update the Newton's guess
            u1_new = u1_new + u1_increment;
            
            
            % Evaluate the norm of the Newton's increment and check if
            % we are close to convergence
            increment_norm = norm(u1_increment(index_u1), 2);
            
            if (increment_norm < tolNewtonsMethod_u)
                flag_isCloseToConvergence = 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over Newton's method
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   2. Find the displacement field (u2)
        %   
        %   Begin: Loop over Newton's method (k)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % By definition, the initial Newton's guess is the solution vector
        % from the previous alternation step
        u2_new = u2;
        
        % Initialize the Newton's increment
        u2_increment = vector_of_zeros;
        
        % External force vector due to body force and tractions
        f2_external = vector_of_zeros;
        f2_external(BCF_array2(:, 1)) = BCF_array2(:, 2);
        
        
        % Flag that indicates whether we are close to convergence
        flag_isCloseToConvergence = 0;
        
        % Norm of the residual vector
        residual_norm = 0;
        
        
        for k = 0 : maxNewtonsMethod
            % Initialize the internal force
            f2_internal = vector_of_zeros;
            
            % Initialize the indices of the matrix entries
            index_for_K22 = (1 : numMatrixEntriesPerElement)';
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over elements (e = e1)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Counters for the quadrature point
            q = 1;
            
            
            for e = 1 : numElements1
                % Find the Bezier extraction matrix
                bezierExtractions1_e = bezierExtractions1(:, :, e);
                
                % Evaluate the map derivative dt/dxi (constant)
                dt_dxi = 1 / elementSizes1(e);
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the parametric domain
                %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                % Find the positions of the nodes
                nodes_e = nodes(IEN_array(:, e), :)';
                
                % Evaluate the B-splines
                basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
                basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the map derivatives (x = x1)
                %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                dx_dxi = nodes_e * basis_der1;
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over quadrature points (q, q_e)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Initialize the element tangent matrix and the element
                % internal force vector
                K22_e = zeros(numDOFsPerElement);
                f2_e  = zeros(numDOFsPerElement, 1);
                
                % Find the equation indices for the element
                index_equation1 = LM_array1(:, e);
                index_equation2 = LM_array2(:, e);
                
                % Get the current Newton's guess for the solution vector
                u1_e = u1(index_equation1)';
                u2_e = u2_new(index_equation2)';
                
                
                for q_e = 1 : numQuadraturePointsPerElement
                    %------------------------------------------------------
                    %  Form the Jacobian matrix
                    %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                    %------------------------------------------------------
                    JacobianMatrix = dx_dxi(:, q_e);
                    
                    % Evaluate the Jacobian
                    Jacobian = JacobianMatrix / dt_dxi;
                    
                    if (Jacobian <= 0)
                        fprintf('\n');
                        fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                        
                        quit;
                    end
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the basis functions in the physical domain
                    %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                    %------------------------------------------------------
                    JacobianMatrix_inv = 1 / JacobianMatrix;
                    
                    basis_physical_der0 =                      basis_der0(:, q_e);
                    basis_physical_der1 = JacobianMatrix_inv * basis_der1(:, q_e);
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the displacement field
                    %------------------------------------------------------
                    u_der1 = u1_e * basis_physical_der1;
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the phase field
                    %------------------------------------------------------
                    c_der0 = u2_e * basis_physical_der0;
                    c_der1 = u2_e * basis_physical_der1;
                    
                    
                    %------------------------------------------------------
                    %  Form the element matrix K22
                    %------------------------------------------------------
                    % Some useful constants for augmented Lagrangian method
                    constant0 = constant_stiffness * u_der1^2;
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscripts denote the derivatives of N_A and N_B.
                    C00 = constant0 + constant_G_c;
                    C11 = constant_G_c * (4 * constant_ell_0_sq);
                    
                    
                    % Additional terms from the augmented Lagrangian method
                    C00 = C00 + lambda0(q) * constant0 + LMM_beta * (3/2 * constant0^2 * c_der0^2);
                    
                    
                    K22_e = K22_e + w(q_e) * Jacobian * (C00 * (basis_physical_der0 * basis_physical_der0') ...
                                                       + C11 * (basis_physical_der1 * basis_physical_der1'));
                    
                    
                    %------------------------------------------------------
                    %  Form the element RHS vector f2
                    %------------------------------------------------------
                    % Some useful constants for augmented Lagrangian method
                    constant0 = constant_stiffness * c_der0 * u_der1^2;
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscript denotes the derivative of N_A.
                    C0 = constant0 - constant_G_c * (1 - c_der0);
                    C1 = constant_G_c * (4 * constant_ell_0_sq * c_der1);
                    
                    
                    % Additional terms from the augmented Lagrangian method
                    C0 = C0 + lambda0(q) * constant0 + LMM_beta * (1/2 * constant0^2 * c_der0);
                    
                    
                    f2_e = f2_e + w(q_e) * Jacobian * (C0 * basis_physical_der0 ...
                                                     + C1 * basis_physical_der1);
                    
                    
                    % Increment the counter for quadrature point
                    q = q + 1;
                end
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Loop over quadrature points
                % ---------------------------------------------------------
                %----------------------------------------------------------
                
                
                %----------------------------------------------------------
                %  Global assembly
                %----------------------------------------------------------
                % Add the element matrix K22
                rows_for_K22(index_for_K22)    = kron(vector_of_ones, index_equation2);
                columns_for_K22(index_for_K22) = kron(index_equation2, vector_of_ones);
                values_for_K22(index_for_K22)  = reshape(K22_e, numMatrixEntriesPerElement, 1);
                index_for_K22 = index_for_K22 + numMatrixEntriesPerElement;
                
                
                % Add the element RHS vector f2
                f2_internal(index_equation2) = f2_internal(index_equation2) + f2_e;
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over elements
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Check for convergence
            %  
            %  We accept that the Newton's guess corresponds to an
            %  equilibrium state via a two-step check.
            %--------------------------------------------------------------
            if (k > 0)
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
                    fprintf('  Equation 2 (phase field equation) has been solved.\n\n');
                    
                    % Update the phase field
                    u2 = u2_new;
                    
                    break;
                    
                elseif (k == maxNewtonsMethod)
                    fprintf('\n');
                    fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                    fprintf('  Equation 2 (phase field equation) has been solved.\n\n');
                    
                    % Update the phase field
                    u2 = u2_new;
                    
                    break;
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Solve for the Newton's increment
            %--------------------------------------------------------------
            % Assemble the tangent matrix
            K22 = sparse(rows_for_K22, columns_for_K22, values_for_K22, numDOFs, numDOFs);
            
            
            % Save the norm of the residual vector at the previous iteration
            residual_norm0 = residual_norm;
            
            % Form the residual vector (note the negative sign)
            residual = f2_external - f2_internal;
            
            % Evaluate the norm of the residual vector at the current iteration
            residual_norm = norm(residual(index_u2), 2);
            
            
            % Apply a row preconditioner
            precond = spdiags(1./max(K22, [], 2), 0, numDOFs, numDOFs);
            K22 = precond * K22;
            residual = precond * residual;
            
            
            % Upon initialization, we told the solution vectors to satisfy
            % the displacements prescribed on the boundary. Hence, we tell
            % the Newton's increment to satisfy zero displacements on the
            % boundary.
            u2_increment(index_u2) = K22(index_u2, index_u2) \ residual(index_u2);
            
            
            % Update the Newton's guess
            u2_new = u2_new + u2_increment;
            
            
            % Evaluate the norm of the Newton's increment and check if
            % we are close to convergence
            increment_norm = norm(u2_increment(index_u2), 2);
            
            if (increment_norm < tolNewtonsMethod_c)
                flag_isCloseToConvergence = 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over Newton's method
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   3. Update the Lagrange multiplier fields and penalty parameter
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counters for the quadrature point
        q = 1;
        
        
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
            %  Loop over quadrature points
            %--------------------------------------------------------------
            % Find the coefficients of the displacement and phase fields
            u1_e = u1(LM_array1(:, e))';
            u2_e = u2(LM_array2(:, e))';
            
            
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
                
                
                %----------------------------------------------------------
                %  Update the Lagrange multiplier fields
                %----------------------------------------------------------
                % Some useful constants for augmented Lagrangian method
                constant0 = 1/2 * constant_stiffness * u_der1^2;
                
                lambda0(q) = lambda0(q) + LMM_beta * constant0 * c_der0^2;
                
                
                % Increment the counter for quadrature point
                q = q + 1;
            end
        end
        
        
        % Update the penalty parameter
        LMM_beta = 10 * LMM_beta;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   4. Save the results
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        save(sprintf('%sfile_results_alternation%d', path_to_results_directory, alternation), 'u1', 'u2', 'lambda0', '-v6');
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over alterations
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end