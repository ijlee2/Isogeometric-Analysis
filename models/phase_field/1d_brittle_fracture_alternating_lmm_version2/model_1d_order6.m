%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 1D cracked bar in tension, so that
%  we can validate the 6th-order phase field theory for brittle fracture.
%  
%  We alternatingly solve the momentum equation (direct method) and the
%  phase field equation (with Newton's method) to get the displacement
%  and phase fields. Lagrange multipliers are added to the phase field
%  equation, so that we can try to enforce the continuity requirement of
%  the phase field at the crack.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_order6.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      restartTime is the alternation step at which we start the simulation
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d_order6(path_to_assembly_directory, path_to_results_directory, restartTime)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  6th-order phase field theory in 1D.\n\n');
    
    
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
         'maxAlternations'    , ...
         'maxNewtonsMethod'   , ...
         'tolNewtonsMethod');
    
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
         'LM_array1'      , 'LM_array2'      , ...
         'BCU_array1'     , 'BCU_array2'     , ...
         'BCF_array1'     , 'BCF_array2'     , ...
         'numUnknownDOFs1', 'numUnknownDOFs2', ...
         'index_u1'       , 'index_u2'       , ...
         'index_f1'       , 'index_f2');
    
    
    
    %----------------------------------------------------------------------
    %  Initialize the fields
    %----------------------------------------------------------------------
    if (restartTime == 0)
        initialize(6, path_to_assembly_directory, path_to_results_directory);
        load(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda');
        
    else
        load(sprintf('%sfile_results_alternation%d', path_to_results_directory, restartTime), 'u1', 'u2', 'lambda');
        
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
    Bernstein1_der2 = zeros(constant_p1p1, numQuadraturePoints);
    Bernstein1_der3 = zeros(constant_p1p1, numQuadraturePoints);
    
    for q_e = 1 : numQuadraturePoints
        temp = eval_1d_bernstein_der(z1(q_e), p1);
        
        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
        Bernstein1_der2(:, q_e) = temp(:, 3);
        Bernstein1_der3(:, q_e) = temp(:, 4);
    end
    
    clear z1 temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over alterations
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_stiffness = material_E * material_A;
    constant_G_c       = material_G_c * material_A / (2 * material_ell_0);
    constant_ell_0_sq  = material_ell_0^2;
    
    
    % Total number of DOFs for the momentum and phase field equations
    numConstraints = 1;
    numDOFs_total1 = numDOFs;
    numDOFs_total2 = numDOFs + numConstraints;
    
    % Some useful constant for global assembly
    constant_p1p1_sq               = constant_p1p1^2;
    constant_p1p1_x_numConstraints = constant_p1p1 * numConstraints;
    
    
    % Initialization for solution and RHS vectors
    globalVector1  = zeros(numDOFs_total1, 1);
    globalVector2  = zeros(numDOFs_total2, 1);
    
    % Initialization for element matrices K11, K12, K22
    elementMatrix1 = zeros(constant_p1p1);
    
    % Initialization for element matrices K13, K23
    elementMatrix2 = zeros(constant_p1p1, numConstraints);
    
    % Initialization for element vectors f1, f2
    elementVector1 = zeros(constant_p1p1, 1);
    
    % Initialization for element vector f3
    elementVector2 = zeros(numConstraints, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones1 = ones(constant_p1p1, 1);
    vector_of_ones2 = ones(numConstraints, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix1 = (1 : constant_p1p1_sq              )';
    indices_for_elementMatrix2 = (1 : constant_p1p1_x_numConstraints)';
    
    
    % Number of matrix entries that we compute for the momentum and
    % phase field equations
    numMatrixEntriesPerElement1 = constant_p1p1_sq;
    numMatrixEntriesPerElement2 = constant_p1p1_sq + 2 * constant_p1p1_x_numConstraints;
    
    % Initialize row, column, value arrays for the stiffness matrix K11
    % and the tangent matrix [K22, K23; K32, K33]
    temp1 = zeros(numMatrixEntriesPerElement1 * numElements1, 1);
    temp2 = zeros(numMatrixEntriesPerElement2 * numElements1, 1);
    
    rows_for_K11    = temp1;
    columns_for_K11 = temp1;
    values_for_K11  = temp1;
    
    rows_for_K22    = temp2;
    columns_for_K22 = temp2;
    values_for_K22  = temp2;
    
    clear temp1 temp2;
    
    
    for alternation = (restartTime + 1) : maxAlternations
        fprintf('\n');
        fprintf('- Alternation index = %d\n', alternation);
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   1. Find the displacement field (u1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Initialize the RHS vector and prescribe the known forces
        f1 = globalVector1;
        f1(BCF_array1(:, 1)) = BCF_array1(:, 2);
        
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
            basis_der0 =             bezierExtractions1_e  * Bernstein1_der0;
            basis_der1 = (dt_dxi   * bezierExtractions1_e) * Bernstein1_der1;
            basis_der2 = (dt_dxi^2 * bezierExtractions1_e) * Bernstein1_der2;
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx_dxi = nodes_e * basis_der1;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over quadrature points (q_e)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Initialize the element stiffness matrix and the element
            % internal force vector
            K11_e = elementMatrix1;
%           f1_e  = elementVector1;
            
            % Find the equation indices for the element
            index_equation1 = LM_array1(:, e);
            index_equation2 = LM_array2(:, e);
            
            % Find the coefficients of the phase field
            u2_e = u2(index_equation2)';
            u3_e = lambda';
            
            
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
                
                basis_physical_der0 =                        basis_der0(:, q_e);
                basis_physical_der1 = JacobianMatrix_inv   * basis_der1(:, q_e);
                basis_physical_der2 = JacobianMatrix_inv^2 * basis_der2(:, q_e);
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                c_der0 = u2_e * basis_physical_der0;
                c_der2 = u2_e * basis_physical_der2;
                
                
                %----------------------------------------------------------
                %  Form the element matrix K11
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
%                C11 = constant_stiffness * c_der0^2;
                C11 = constant_stiffness * (c_der0^2 + u3_e(1) * (material_ell_0^2 * c_der2)^2);
                
                K11_e = K11_e + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1');
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f1
                %----------------------------------------------------------
%               f1_e = f1_e + (w(q_e) * Jacobian * 0) * basis_physical_der0;
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Global assembly
            %--------------------------------------------------------------
            % Add the element matrix K11
            index = lastIndex_for_K + indices_for_elementMatrix1;
            
            rows_for_K11(index)    = kron(vector_of_ones1, index_equation1);
            columns_for_K11(index) = kron(index_equation1, vector_of_ones1);
            values_for_K11(index)  = reshape(K11_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element RHS vector f1
%           f1(index_equation1) = f1(index_equation1) + f1_e;
            
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Solve for the displacement field
        %------------------------------------------------------------------
        % Assemble the stiffness matrix
        K11 = sparse(rows_for_K11, columns_for_K11, values_for_K11, numDOFs_total1, numDOFs_total1);
        
        
        % Apply a row preconditioner
        precond = spdiags(1./max(K11, [], 2), 0, numDOFs_total1, numDOFs_total1);
        K11 = precond * K11;
        f1  = precond * f1;
        
        
        % Prescribe the known coefficients and solve for the unknown
        u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
        u1(index_u1) = K11(index_u1, index_u1) \ (f1(index_u1) - K11(index_u1, index_f1) * u1(index_f1));
        
        
        fprintf('\n');
        fprintf('  Equation 1 (momentum equation) has been solved.\n\n');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   2. Find the phase field (u2)
        %   
        %   Begin: Loop over Newton's method (k)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Set the Lagrange multipliers to 0
        lambda = zeros(numConstraints, 1);
        
        % By definition, the initial Newton's guess is the solution vector
        % from the previous time step
        u_new = [u2; lambda];
        
        % External force vector due to body force and tractions
        f_external = globalVector2;
        f_external(BCF_array2(:, 1)) = BCF_array2(:, 2);
        
        
        % Initialize the Newton's increment
        u_increment = globalVector2;
        
        % Indices of the unknown coefficients
        index_equation3 = numDOFs + (1 : numConstraints)';
        index_u = [index_u2; index_equation3];
        
        
        % Flag that indicates whether we are close to convergence
        flag_isCloseToConvergence = 0;
        
        % Norm of the residual vector
        residual_norm = 0;
        
        
        for k = 0 : maxNewtonsMethod
            % Initialize the internal force
            f_internal = globalVector2;
            
            % Index that was last used to set the entry of K
            lastIndex_for_K = 0;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over elements (e = e1)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
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
                basis_der0 =             bezierExtractions1_e  * Bernstein1_der0;
                basis_der1 = (dt_dxi   * bezierExtractions1_e) * Bernstein1_der1;
                basis_der2 = (dt_dxi^2 * bezierExtractions1_e) * Bernstein1_der2;
                basis_der3 = (dt_dxi^3 * bezierExtractions1_e) * Bernstein1_der3;
                
                
                %----------------------------------------------------------
                %  Evaluate the map derivatives (x = x1)
                %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                dx_dxi = nodes_e * basis_der1;
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over quadrature points (q_e)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Initialize the element tangent matrices and the element
                % internal force vectors
                K22_e = elementMatrix1;
                K23_e = elementMatrix2;
                
                f2_e  = elementVector1;
                f3_e  = elementVector2;
                
                % Find the equation indices for the element
                index_equation1 = LM_array1(:, e);
                index_equation2 = LM_array2(:, e);
                
                % Find the coefficients of the displacement field
                u1_e = u1(index_equation1)';
                
                % Get the current Newton's guess for the phase field and
                % Lagrange multipliers
                u2_e = u_new(index_equation2)';
                u3_e = u_new(index_equation3)';
                
                
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
                    
                    basis_physical_der0 =                        basis_der0(:, q_e);
                    basis_physical_der1 = JacobianMatrix_inv   * basis_der1(:, q_e);
                    basis_physical_der2 = JacobianMatrix_inv^2 * basis_der2(:, q_e);
                    basis_physical_der3 = JacobianMatrix_inv^3 * basis_der3(:, q_e);
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the displacement field
                    %------------------------------------------------------
                    u_der1 = u1_e * basis_physical_der1;
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the phase field
                    %------------------------------------------------------
                    c_der0 = u2_e * basis_physical_der0;
                    c_der1 = u2_e * basis_physical_der1;
                    c_der2 = u2_e * basis_physical_der2;
                    c_der3 = u2_e * basis_physical_der3;
                    
                    
                    %------------------------------------------------------
                    %  Form the element matrices K22, K23
                    %------------------------------------------------------
                    % Some useful constants for Lagrange multiplier method
                    constant0 = constant_stiffness                       * u_der1^2;
                    constant2 = constant_stiffness * constant_ell_0_sq^2 * u_der1^2;
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscripts denote the derivatives of N_A and N_B.
                    C00 = constant0 + constant_G_c;
                    C11 =             constant_G_c * (4/3    * constant_ell_0_sq  );
                    C22 =             constant_G_c * (16/27  * constant_ell_0_sq^2);
                    C33 =             constant_G_c * (64/729 * constant_ell_0_sq^3);
                    
                    
                    % Additional terms from Lagrange multiplier method
                    C22 = C22 + u3_e(1) * constant2;
                    
                    
                    K22_e = K22_e + (w(q_e) * Jacobian * C00) * (basis_physical_der0 * basis_physical_der0') ...
                                  + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1') ...
                                  + (w(q_e) * Jacobian * C22) * (basis_physical_der2 * basis_physical_der2') ...
                                  + (w(q_e) * Jacobian * C33) * (basis_physical_der3 * basis_physical_der3');
                    
                    K23_e(:, 1) = K23_e(:, 1) + (w(q_e) * Jacobian * constant2 * c_der2) * basis_physical_der2;
                    
                    
                    %------------------------------------------------------
                    %  Form the element RHS vectors f2, f3
                    %------------------------------------------------------
                    % Some useful constants for Lagrange multiplier method
                    constant0 = constant_stiffness                       * c_der0 * u_der1^2;
                    constant2 = constant_stiffness * constant_ell_0_sq^2 * c_der2 * u_der1^2;
                    
                    vector_of_constants = [0.5 * constant2 * c_der2];
                    
                    
                    % Coefficients for the (regular) phase field theory.
                    % The subscript denotes the derivative of N_A.
                    C0 = constant0 - constant_G_c * (1 - c_der0);
                    C1 =             constant_G_c * (4/3    * constant_ell_0_sq   * c_der1);
                    C2 =             constant_G_c * (16/27  * constant_ell_0_sq^2 * c_der2);
                    C3 =             constant_G_c * (64/729 * constant_ell_0_sq^3 * c_der3);
                    
                    
                    % Additional terms from Lagrange multiplier method
                    C2 = C2 + u3_e(1) * constant2;
                    
                    
                    f2_e = f2_e + (w(q_e) * Jacobian * C0) * basis_physical_der0 ...
                                + (w(q_e) * Jacobian * C1) * basis_physical_der1 ...
                                + (w(q_e) * Jacobian * C2) * basis_physical_der2 ...
                                + (w(q_e) * Jacobian * C3) * basis_physical_der3;
                    
                    f3_e = f3_e + (w(q_e) * Jacobian) * vector_of_constants;
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
                index = lastIndex_for_K + indices_for_elementMatrix1;
                
                rows_for_K22(index)    = kron(vector_of_ones1, index_equation2);
                columns_for_K22(index) = kron(index_equation2, vector_of_ones1);
                values_for_K22(index)  = reshape(K22_e, constant_p1p1_sq, 1);
                
                lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
                
                
                % Add the element matrix K23
                index = lastIndex_for_K + indices_for_elementMatrix2;
                
                rows_for_K22(index)    = kron(vector_of_ones2, index_equation2);
                columns_for_K22(index) = kron(index_equation3, vector_of_ones1);
                values_for_K22(index)  = reshape(K23_e, constant_p1p1_x_numConstraints, 1);
                
                lastIndex_for_K = lastIndex_for_K + constant_p1p1_x_numConstraints;
                
                
                % Add the element matrix K32
                index = lastIndex_for_K + indices_for_elementMatrix2;
                
                rows_for_K22(index)    = kron(vector_of_ones1, index_equation3);
                columns_for_K22(index) = kron(index_equation2, vector_of_ones2);
                values_for_K22(index)  = reshape(K23_e', constant_p1p1_x_numConstraints, 1);
                
                lastIndex_for_K = lastIndex_for_K + constant_p1p1_x_numConstraints;
                
                
                % Add the element matrix K33
                % Do nothing (K33 = 0)
                
                
                % Add the element RHS vector f2
                f_internal(index_equation2) = f_internal(index_equation2) + f2_e;
                
                
                % Add the element RHS vector f3
                f_internal(index_equation3) = f_internal(index_equation3) + f3_e;
                
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over elements
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Check for convergence
            %  
            %  We accept that the Newton's guess corresponds to an equilibrium
            %  state via a two-step check.
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
                    
                    % Update the phase field and Lagrange multipliers
                    u2     = u_new((1 : numDOFs)');
                    lambda = u_new(index_equation3);
                    
                    break;
                    
                elseif (k == maxNewtonsMethod)
                    fprintf('\n');
                    fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                    
                    % Update the phase field and Lagrange multipliers
                    u2     = u_new((1 : numDOFs)');
                    lambda = u_new(index_equation3);
                    
                    break;
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Solve for the Newton's increment
            %--------------------------------------------------------------
            % Assemble the tangent matrix
            K22 = sparse(rows_for_K22, columns_for_K22, values_for_K22, numDOFs_total2, numDOFs_total2);
            
            
            % Save the norm of the residual vector at the previous iteration
            residual_norm0 = residual_norm;
            
            % Form the residual vector (note the negative sign)
            residual = f_external - f_internal;
            
            % Evaluate the norm of the residual vector at the current iteration
            residual_norm = norm(residual(index_u2), 2);
            
            
            % Apply a row preconditioner
            precond = spdiags(1./max(K22, [], 2), 0, numDOFs_total2, numDOFs_total2);
            K22      = precond * K22;
            residual = precond * residual;
            
            
            % Upon initialization, we told the solution vectors to satisfy the
            % displacements prescribed on the boundary. Hence, we tell the
            % Newton's increment to satisfy zero displacements on the boundary.
            u_increment(index_u) = K22(index_u, index_u) \ residual(index_u);
            
            
            % Update the Newton's guess
            u_new = u_new + u_increment;
            
            
            % Evaluate the norm of the Newton's increment and check if
            % we are close to convergence
            increment_norm = norm(u_increment(index_u2), 2);
            
            if (increment_norm < tolNewtonsMethod)
                flag_isCloseToConvergence = 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over Newton's method
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        fprintf('\n');
        fprintf('  Equation 2 (phase field equation) has been solved.\n\n');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   3. Save the results
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        save(sprintf('%sfile_results_alternation%d', path_to_results_directory, alternation), 'u1', 'u2', 'lambda', '-v6');
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
