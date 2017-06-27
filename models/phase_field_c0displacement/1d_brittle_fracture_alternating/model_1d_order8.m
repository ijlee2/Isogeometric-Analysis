%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 1D cracked bar in tension, so that
%  we can validate the 8th-order phase field theory for brittle fracture.
%  
%  We alternatingly solve the momentum and phase field equations (direct
%  method) to get the displacement and phase fields.
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
%      restartTime is the alternation step at which we start the simulation
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d_order8(path_to_assembly_directory, path_to_results_directory, restartTime)
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
         'material_A'         , ...
         'material_G_c'       , ...
         'material_ell_0'     , ...
         ...
         'maxAlternations');
    
    
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
        initialize(8, path_to_assembly_directory, path_to_results_directory);
        load(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', 'u2');
        
    else
        load(sprintf('%sfile_results_alternation%d', path_to_results_directory, restartTime), 'u1', 'u2');
        
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
    Bernstein1_der4 = zeros(constant_p1p1, numQuadraturePoints);
    
    for q_e = 1 : numQuadraturePoints
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
    constant_stiffness = material_E * material_A;
    constant_G_c       = material_G_c * material_A / (2 * material_ell_0);
    constant_ell_0_sq  = material_ell_0^2;
    
    % Some useful constant for global assembly
    constant_p1p1_sq = constant_p1p1^2;
    
    
    % Initialization for solution and RHS vectors
    globalVector1 = zeros(numDOFs(1), 1);
    globalVector2 = zeros(numDOFs(2), 1);
    
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
    numMatrixEntriesPerElement1 = constant_p1p1_sq;
    numMatrixEntriesPerElement2 = constant_p1p1_sq;
    
    % Initialize row, column, value arrays for the stiffness matrices
    % K11 and K22
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
%           basis1_der0 =           bezierExtractions11_e  * Bernstein1_der0;
            basis1_der1 = (dt_dxi * bezierExtractions11_e) * Bernstein1_der1;
            
            % Evaluate the B-splines for the phase field
            basis2_der0 =           bezierExtractions21_e  * Bernstein1_der0;
            basis2_der1 = (dt_dxi * bezierExtractions21_e) * Bernstein1_der1;
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx1_dxi = nodes1_e * basis1_der1;
            dx2_dxi = nodes2_e * basis2_der1;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over quadrature points (q_e)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Initialize the element stiffness matrix and the element
            % internal force vector
            K11_e = elementMatrix;
%           f1_e  = elementVector;
            
            % Find the equation indices for the element
            index_equation1 = LM_array1(:, e);
            index_equation2 = LM_array2(:, e);
            
            % Find the coefficients of the phase field
            u2_e = u2(index_equation2)';
            
            
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
%               JacobianMatrix2_inv = 1 / JacobianMatrix2;
                
                % For the displacement field
%               basis1_physical_der0 =                       basis1_der0(:, q_e);
                basis1_physical_der1 = JacobianMatrix1_inv * basis1_der1(:, q_e);
                
                % For the phase field
                basis2_physical_der0 =                       basis2_der0(:, q_e);
%               basis2_physical_der1 = JacobianMatrix2_inv * basis2_der1(:, q_e);
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                c_der0 = u2_e * basis2_physical_der0;
                
                
                %----------------------------------------------------------
                %  Form the element matrix K11
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
                C11 = constant_stiffness * c_der0^2;
                
                K11_e = K11_e + (w(q_e) * Jacobian1 * C11) * (basis1_physical_der1 * basis1_physical_der1');
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f1
                %----------------------------------------------------------
%               f1_e = f1_e + (w(q_e) * Jacobian1 * 0) * basis1_physical_der0;
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
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K11(index)    = kron(vector_of_ones, index_equation1);
            columns_for_K11(index) = kron(index_equation1, vector_of_ones);
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
        K11 = sparse(rows_for_K11, columns_for_K11, values_for_K11, numDOFs(1), numDOFs(1));
        
        
        % Apply a row preconditioner
        precond = spdiags(1./max(K11, [], 2), 0, numDOFs(1), numDOFs(1));
        K11 = precond * K11;
        f1  = precond * f1;
        
        
        % Prescribe the known coefficients and solve for the unknown
        u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
        u1(index_u1) = K11(index_u1, index_u1) \ (f1(index_u1) - K11(index_u1, index_f1) * u1(index_f1));
        
        
        fprintf('\n');
        fprintf('  Equation 1 (momentum equation) has been solved.\n');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   2. Find the phase field (u2)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Initialize the RHS vector and prescribe the known forces
        f2 = globalVector2;
        f2(BCF_array2(:, 1)) = BCF_array2(:, 2);
        
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
            basis2_der4 = (dt_dxi^4 * bezierExtractions21_e) * Bernstein1_der4;
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx1_dxi = nodes1_e * basis1_der1;
            dx2_dxi = nodes2_e * basis2_der1;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over quadrature points (q_e)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % Initialize the element stiffness matrix and the element
            % internal force vector
            K22_e = elementMatrix;
            f2_e  = elementVector;
            
            % Find the equation indices for the element
            index_equation1 = LM_array1(:, e);
            index_equation2 = LM_array2(:, e);
            
            % Find the coefficients of the displacement field
            u1_e = u1(index_equation1)';
            
            
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
                basis2_physical_der4 = JacobianMatrix2_inv^4 * basis2_der4(:, q_e);
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u_der1 = u1_e * basis1_physical_der1;
                
                
                %----------------------------------------------------------
                %  Form the element matrix K22
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscripts denote the derivatives of N_A and N_B.
                C00 = constant_stiffness * u_der1^2 + constant_G_c;
                C11 =                                 constant_G_c * (        constant_ell_0_sq  );
                C22 =                                 constant_G_c * (3/8   * constant_ell_0_sq^2);
                C33 =                                 constant_G_c * (1/16  * constant_ell_0_sq^3);
                C44 =                                 constant_G_c * (1/256 * constant_ell_0_sq^4);
                
                K22_e = K22_e + (w(q_e) * Jacobian2 * C00) * (basis2_physical_der0 * basis2_physical_der0') ...
                              + (w(q_e) * Jacobian2 * C11) * (basis2_physical_der1 * basis2_physical_der1') ...
                              + (w(q_e) * Jacobian2 * C22) * (basis2_physical_der2 * basis2_physical_der2') ...
                              + (w(q_e) * Jacobian2 * C33) * (basis2_physical_der3 * basis2_physical_der3') ...
                              + (w(q_e) * Jacobian2 * C44) * (basis2_physical_der4 * basis2_physical_der4');
                
                
                %----------------------------------------------------------
                %  Form the element RHS vector f2
                %----------------------------------------------------------
                % Coefficients for the (regular) phase field theory.
                % The subscript denotes the derivative of N_A.
                C0 = constant_G_c;
                
                f2_e = f2_e + (w(q_e) * Jacobian2 * C0) * basis2_physical_der0;
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over quadrature points
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Global assembly
            %--------------------------------------------------------------
            % Add the element matrix K22
            index = lastIndex_for_K + indices_for_elementMatrix;
            
            rows_for_K22(index)    = kron(vector_of_ones, index_equation2);
            columns_for_K22(index) = kron(index_equation2, vector_of_ones);
            values_for_K22(index)  = reshape(K22_e, constant_p1p1_sq, 1);
            
            lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
            
            
            % Add the element RHS vector f2
            f2(index_equation2) = f2(index_equation2) + f2_e;
            
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Solve for the phase field
        %------------------------------------------------------------------
        % Assemble the stiffness matrix
        K22 = sparse(rows_for_K22, columns_for_K22, values_for_K22, numDOFs(2), numDOFs(2));
        
        
        % Apply a row preconditioner
        precond = spdiags(1./max(K22, [], 2), 0, numDOFs(2), numDOFs(2));
        K22 = precond * K22;
        f2  = precond * f2;
        
        
        % Prescribe the known coefficients and solve for the unknown
        u2(BCU_array2(:, 1)) = BCU_array2(:, 2);
        u2(index_u2) = K22(index_u2, index_u2) \ (f2(index_u2) - K22(index_u2, index_f2) * u2(index_f2));
        
        
        fprintf('\n');
        fprintf('  Equation 2 (phase field equation) has been solved.\n\n');
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   3. Save the results
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        save(sprintf('%sfile_results_alternation%d', path_to_results_directory, alternation), 'u1', 'u2', '-v6');
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
