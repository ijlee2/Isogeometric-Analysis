%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 1D cracked bar in tension, so that
%  we can validate the 6th-order phase field theory for brittle fracture.
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
%      ./matbg.sh driver_order6.m output
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
    
    % Things specific to the displacement field
    nodes1              = nodes;
    bezierExtractions11 = bezierExtractions1;
    IEN_array1          = IEN_array;
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 2), ...
         'nodes'             , ...
         'bezierExtractions1', ...
         'IEN_array');
    
    % Things specific to the phase field
    nodes2              = nodes;
    bezierExtractions21 = bezierExtractions1;
    IEN_array2          = IEN_array;
    
    clear nodes bezierExtractions1 IEN_array;
    
    
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
        load(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2');
        
    else
        load(sprintf('%sfile_results_iteration%d', path_to_results_directory, restartTime), 'u1', 'u2');
        
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
    %   Begin: Loop over Newton's method (k)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_stiffness = material_E * material_A;
    constant_G_c       = material_G_c * material_A / (2 * material_ell_0);
    constant_ell_0_sq  = material_ell_0^2;
    
    
    % Total number of DOFs due to the displacement and phase fields
    numDOFs_total = sum(numDOFs);
    
    % Some useful constant for global assembly
    constant_p1p1_sq = constant_p1p1^2;
    
    
    % Initialization for solution and RHS vectors
    globalVector   = zeros(numDOFs_total, 1);
    
    % Initialization for element matrices K11, K12, K22
    elementMatrix = zeros(constant_p1p1);
    
    % Initialization for element vectors f1, f2
    elementVector = zeros(constant_p1p1, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(constant_p1p1, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : constant_p1p1_sq)';
    
    
    % Number of matrix entries that we compute for the tangent matrix Ku3_e(1)
    numMatrixEntriesPerElement = 4 * constant_p1p1_sq;
    
    % Initialize row, column, value arrays for the tangent matrix K
    temp = zeros(numMatrixEntriesPerElement * numElements1, 1);
    
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
    
    clear temp;
    
    
    % By definition, the initial Newton's guess is the solution vector
    % from the previous time step
    u_new = [u1; u2];
    
    % External force vector due to increments in body force and tractions
    f_external = globalVector;
    f_external(BCF_array1(:, 1))              = BCF_array1(:, 2);
    f_external(BCF_array2(:, 1) + numDOFs(1)) = BCF_array2(:, 2);
    
    
    % Initialize the Newton's increment
    u_increment = globalVector;
    
    % Indices of the unknown coefficients
    index_u = [index_u1; index_u2 + numDOFs(1)];
    
    
    % Flag that indicates whether we are close to convergence
    flag_isCloseToConvergence = 0;
    
    % Norm of the residual vector
    residual_norm = 0;
    
    
    for k = restartTime : maxNewtonsMethod
        % Initialize the internal force
        f_internal = globalVector;
        
        % Index that was last used to set the entry of K
        lastIndex_for_K = 0;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        for e = 1 : numElements1
            % Find the Bezier extraction matrix
            bezierExtractions11_e = bezierExtractions11(:, :, e);
            bezierExtractions21_e = bezierExtractions21(:, :, e);
            
            % Evaluate the map derivative dt/dxi (constant)
            dt_dxi = 1 / elementSizes1(e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the parametric domain
            %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            % Find the positions of the nodes
            nodes1_e = nodes1(IEN_array1(:, e), :)';
            nodes2_e = nodes2(IEN_array2(:, e), :)';
            
            % Evaluate the B-splines for the displacement field
%           basis1_der0 =             bezierExtractions11_e  * Bernstein1_der0;
            basis1_der1 = (dt_dxi   * bezierExtractions11_e) * Bernstein1_der1;
            
            % Evaluate the B-splines for the phase field
            basis2_der0 =             bezierExtractions21_e  * Bernstein1_der0;
            basis2_der1 = (dt_dxi   * bezierExtractions21_e) * Bernstein1_der1;
            basis2_der2 = (dt_dxi^2 * bezierExtractions21_e) * Bernstein1_der2;
            basis2_der3 = (dt_dxi^3 * bezierExtractions21_e) * Bernstein1_der3;
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx1_dxi = nodes1_e * basis1_der1;
            dx2_dxi = nodes2_e * basis2_der1;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over quadrature points (q, q_e)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Initialize the element tangent matrices and the element
            % internal force vectors
            K11_e = elementMatrix;
            K12_e = elementMatrix;
            K22_e = elementMatrix;
            
            f1_e  = elementVector;
            f2_e  = elementVector;
            
            % Find the equation indices for the element
            index_equation1 = LM_array1(:, e);
            index_equation2 = LM_array2(:, e) + numDOFs(1);
            
            % Get the current Newton's guess for the displacement and
            % phase fields
            u1_e = u_new(index_equation1)';
            u2_e = u_new(index_equation2)';
            
            
            for q_e = 1 : numQuadraturePointsPerElement
                %----------------------------------------------------------
                %  Form the Jacobian matrix
                %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                %----------------------------------------------------------
                JacobianMatrix1 = dx1_dxi(:, q_e);
                JacobianMatrix2 = dx2_dxi(:, q_e);
                
                % Evaluate the Jacobian
                Jacobian1 = JacobianMatrix1 / dt_dxi;
                Jacobian2 = JacobianMatrix2 / dt_dxi;
                
                if (Jacobian1 <= 0 || Jacobian2 <= 0)
                    fprintf('\n');
                    fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                    
                    quit;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the physical domain
                %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                %----------------------------------------------------------
                JacobianMatrix1_inv = 1 / JacobianMatrix1;
                JacobianMatrix2_inv = 1 / JacobianMatrix2;
                
                % For the displacement field
%               basis1_physical_der0 =                         basis1_der0(:, q_e);
                basis1_physical_der1 = JacobianMatrix1_inv   * basis1_der1(:, q_e);
                
                % For the phase field
                basis2_physical_der0 =                         basis2_der0(:, q_e);
                basis2_physical_der1 = JacobianMatrix2_inv   * basis2_der1(:, q_e);
                basis2_physical_der2 = JacobianMatrix2_inv^2 * basis2_der2(:, q_e);
                basis2_physical_der3 = JacobianMatrix2_inv^3 * basis2_der3(:, q_e);
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u_der1 = u1_e * basis1_physical_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                c_der0 = u2_e * basis2_physical_der0;
                c_der1 = u2_e * basis2_physical_der1;
                c_der2 = u2_e * basis2_physical_der2;
                c_der3 = u2_e * basis2_physical_der3;
                
                
                %----------------------------------------------------------
                %  Form the element matrix K11
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
                C11 = constant_stiffness * c_der0^2;
                
                
                K11_e = K11_e + (w(q_e) * Jacobian1 * C11) * (basis1_physical_der1 * basis1_physical_der1');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K12
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
                C10 = 2 * constant_stiffness * c_der0 * u_der1;
                
                
                K12_e = K12_e + (w(q_e) * Jacobian1 * C10) * (basis1_physical_der1 * basis2_physical_der0');
                
                
                %----------------------------------------------------------
                %  Form the element matrix K22
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
                C00 = constant_stiffness * u_der1^2 + constant_G_c                                 ;
                C11 =                                 constant_G_c * (4/3    * constant_ell_0_sq  );
                C22 =                                 constant_G_c * (16/27  * constant_ell_0_sq^2);
                C33 =                                 constant_G_c * (64/729 * constant_ell_0_sq^3);
                
                
                K22_e = K22_e + (w(q_e) * Jacobian2 * C00) * (basis2_physical_der0 * basis2_physical_der0') ...
                              + (w(q_e) * Jacobian2 * C11) * (basis2_physical_der1 * basis2_physical_der1') ...
                              + (w(q_e) * Jacobian2 * C22) * (basis2_physical_der2 * basis2_physical_der2') ...
                              + (w(q_e) * Jacobian2 * C33) * (basis2_physical_der3 * basis2_physical_der3');
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f1
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscript denotes the derivative of N_A.
                C1 = constant_stiffness * c_der0^2 * u_der1;
                
                
                f1_e = f1_e + (w(q_e) * Jacobian1 * C1) * basis1_physical_der1;
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f2
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscript denotes the derivative of N_A.
                C0 = constant_stiffness * c_der0 * u_der1^2 - constant_G_c * (1 - c_der0)                           ;
                C1 =                                          constant_G_c * (4/3    * constant_ell_0_sq   * c_der1);
                C2 =                                          constant_G_c * (16/27  * constant_ell_0_sq^2 * c_der2);
                C3 =                                          constant_G_c * (64/729 * constant_ell_0_sq^3 * c_der3);
                
                
                f2_e = f2_e + (w(q_e) * Jacobian2 * C0) * basis2_physical_der0 ...
                            + (w(q_e) * Jacobian2 * C1) * basis2_physical_der1 ...
                            + (w(q_e) * Jacobian2 * C2) * basis2_physical_der2 ...
                            + (w(q_e) * Jacobian2 * C3) * basis2_physical_der3;
            end
            
            
            % Do we set the off-diagonal blocks to zero?
            K12_e = elementMatrix;
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Global assembly
            %--------------------------------------------------------------
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
            
            
            % Add the element RHS vector f1
            f_internal(index_equation1) = f_internal(index_equation1) + f1_e;
            
            
            % Add the element RHS vector f2
            f_internal(index_equation2) = f_internal(index_equation2) + f2_e;
            
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
            %--------------------------------------------------------------
            %  Save the results for debugging
            %--------------------------------------------------------------
            u1 = u_new(1 : numDOFs(1));
            u2 = u_new((numDOFs(1) + 1) : numDOFs_total);
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, k), 'u1', 'u2', '-v6');
            
            
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
                
                % Update the displacement and phase fields
                u1 = u_new(1 : numDOFs(1));
                u2 = u_new((numDOFs(1) + 1) : numDOFs_total);
                
                break;
                
            elseif (k == maxNewtonsMethod)
                fprintf('\n');
                fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                
                % Update the displacement and phase fields
                u1 = u_new(1 : numDOFs(1));
                u2 = u_new((numDOFs(1) + 1) : numDOFs_total);
                
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
        
        % Form the residual vector (note the negative sign)
        residual = f_external - f_internal;
        
        % Evaluate the norm of the residual vector at the current iteration
        residual_norm = norm(residual(index_u), 2);
        
        
        % Apply a row preconditioner
        precond = spdiags(1./max(K, [], 2), 0, numDOFs_total, numDOFs_total);
        K        = precond * K;
        residual = precond * residual;
        
        
        % Upon initialization, we told the solution vectors to satisfy the
        % displacements prescribed on the boundary. Hence, we tell the
        % Newton's increment to satisfy zero displacements on the boundary.
        u_increment(index_u) = K(index_u, index_u) \ residual(index_u);
        
        
        % Update the Newton's guess
        u_new = u_new + u_increment;
        
        
        % Evaluate the norm of the Newton's increment and check if
        % we are close to convergence
        increment_norm = norm(u_increment(index_u), 2);
        
        if (increment_norm < tolNewtonsMethod)
            flag_isCloseToConvergence = 1;
        end
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
