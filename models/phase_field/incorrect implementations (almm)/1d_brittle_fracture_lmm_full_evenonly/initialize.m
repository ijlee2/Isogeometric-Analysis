%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine initializes the coefficients for the displacement and
%  phase field. We perturb the true solutions by a small amount and take
%  the L2 projection to get the coefficients.
%  
%  
%  Instructions:
%  
%  Type the following in a code,
%  
%      initialize(order, path_to_assembly_directory, path_to_results_directory)
%  
%  where,
%  
%      order is the order of the phase field theory (2, 4, 6, 8)
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function initialize(order, path_to_assembly_directory, path_to_results_directory)
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
         'index_f1', 'index_f2'    , ...
         'u_L');
    
    
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
    
    for j = 1 : numQuadraturePoints
        temp = eval_1d_bernstein_der(z1(j), p1);
        
        Bernstein1_der0(:, j) = temp(:, 1);
        Bernstein1_der1(:, j) = temp(:, 2);
    end
    
    clear z1 temp;
    
    
    %----------------------------------------------------------------------
    %  Initialize the L2 projection matrices and the solution and RHS
    %  vectors
    %----------------------------------------------------------------------
    % Row, column, value arrays for the L2 projection matrices K11 and K22
    temp = zeros(numMatrixEntries, 1);
    
    rows_for_K11    = temp;
    columns_for_K11 = temp;
    values_for_K11  = temp;
    
    rows_for_K22    = temp;
    columns_for_K22 = temp;
    values_for_K22  = temp;
    
    clear temp;
    
    
    % Initialize the solution vectors
    u1 = zeros(numDOFs, 1);
    u2 = zeros(numDOFs, 1);
    
    % Initialize the RHS vectors and prescribe the known coefficients
    f1 = zeros(numDOFs, 1);
    f2 = zeros(numDOFs, 1);
    f1(BCF_array1(:, 1)) = BCF_array1(:, 2);
    f2(BCF_array2(:, 1)) = BCF_array2(:, 2);
    
    
    % Number of matrix entries that we compute per element
    numMatrixEntriesPerElement = numDOFsPerElement^2;
    
    % Initialize the indices of the matrix entries
    index_for_K = (1 : numMatrixEntriesPerElement)';
    
    % Vector of ones for global assembly
    vector_of_ones = ones(numDOFsPerElement, 1);
    
    
    
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
        basis_der0 = bezierExtractions1_e * Bernstein1_der0;
        basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
        
        
        %------------------------------------------------------------------
        %  Evaluate the points on the physical domain
        %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
        %------------------------------------------------------------------
        x = nodes_e * basis_der0;
        
        
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
        % Initialize the element matrices and RHS vectors
        K_e = zeros(numDOFsPerElement);
        f1_e = zeros(numDOFsPerElement, 1);
        f2_e = zeros(numDOFsPerElement, 1);
        
        
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
            basis_physical_der0 = basis_der0(:, q_e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the exact solutions
            %--------------------------------------------------------------
            u_exact = function_u_exact(x(q_e), u_L, material_ell_0, order);
            c_exact = function_c_exact(x(q_e), material_ell_0, order);
            
            
            %--------------------------------------------------------------
            %  Add the element L2 projection matrix
            %--------------------------------------------------------------
            K_e = K_e + w(q_e) * Jacobian * (basis_physical_der0 * basis_physical_der0');
            
            
            %--------------------------------------------------------------
            %  Add the element RHS vectors
            %--------------------------------------------------------------
            f1_e = f1_e + w(q_e) * Jacobian * (u_exact * basis_physical_der0);
            f2_e = f2_e + w(q_e) * Jacobian * (c_exact * basis_physical_der0);
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over quadrature points
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Get the global equation indices
        index_equation1 = LM_array1(:, e);
        index_equation2 = LM_array2(:, e);
        
        
        % Add the element matrices
        rows_for_K11(index_for_K)    = kron(vector_of_ones, index_equation1);
        columns_for_K11(index_for_K) = kron(index_equation1, vector_of_ones);
        values_for_K11(index_for_K)  = reshape(K_e, numMatrixEntriesPerElement, 1);
        
        rows_for_K22(index_for_K)    = kron(vector_of_ones, index_equation2);
        columns_for_K22(index_for_K) = kron(index_equation2, vector_of_ones);
        values_for_K22(index_for_K)  = reshape(K_e, numMatrixEntriesPerElement, 1);
        
        index_for_K = index_for_K + numMatrixEntriesPerElement;
        
        
        % Add the element RHS vectors
        f1(index_equation1) = f1(index_equation1) + f1_e;
        f2(index_equation2) = f2(index_equation2) + f2_e;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Solve for the displacement and phase fields
    %----------------------------------------------------------------------
    % Assemble the L2 projection matrices
    K11 = sparse(rows_for_K11, columns_for_K11, values_for_K11, numDOFs, numDOFs);
    K22 = sparse(rows_for_K22, columns_for_K22, values_for_K22, numDOFs, numDOFs);
    
    
    % Apply a row preconditioner
    precond1 = spdiags(1./max(K11, [], 2), 0, numDOFs, numDOFs);
    precond2 = spdiags(1./max(K22, [], 2), 0, numDOFs, numDOFs);
    K11 = precond1 * K11;
    K22 = precond2 * K22;
    f1 = precond1 * f1;
    f2 = precond2 * f2;
    
    
    % Prescribe the known coefficients
    u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
    u2(BCU_array2(:, 1)) = BCU_array2(:, 2);
    
    
    % Solve for the unknown coefficients
    u1(index_u1) = K11(index_u1, index_u1) \ (f1(index_u1) - K11(index_u1, index_f1) * u1(index_f1));
    u2(index_u2) = K22(index_u2, index_u2) \ (f2(index_u2) - K22(index_u2, index_f2) * u2(index_f2));
    
    
    %----------------------------------------------------------------------
    %  Save the fields
    %----------------------------------------------------------------------
    temp = zeros(numQuadraturePointsPerElement * numElements1, 1);
    
    switch order
        case 2
            % Initialize the Lagrange multiplier field
            lambda0 = temp;
            
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda0', '-v6');
            
        case 4
            % Initialize the Lagrange multiplier fields
            lambda0 = temp;
            
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda0', '-v6');
            
        case 6
            % Initialize the Lagrange multiplier fields
            lambda0 = temp;
            lambda2 = temp;
            
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda0', 'lambda2', '-v6');
            
        case 8
            % Initialize the Lagrange multiplier fields
            lambda0 = temp;
            lambda2 = temp;
            
            save(sprintf('%sfile_results_iteration%d', path_to_results_directory, 0), 'u1', 'u2', 'lambda0', 'lambda2', '-v6');
            
    end
end


function u = function_u_exact(x, u_L, epsilon, order)
    switch order
        case 2
            if (x < -epsilon)
                u = 0;
            elseif (x < epsilon)
                u = u_L/(2*epsilon) * (x + epsilon);
            else
                u = u_L;
            end
            
        case 4
            if (x < -epsilon)
                u = 0;
            elseif (x < 0)
                u = u_L/(2*epsilon^2) * (x^2 + 2*epsilon*x + epsilon^2);
            elseif (x < epsilon)
                u = u_L/(2*epsilon^2) * (-x^2 + 2*epsilon*x + epsilon^2);
            else
                u = u_L;
            end
            
        case 6
            if (x < -epsilon)
                u = 0;
            elseif (x < -epsilon/3)
                u = u_L/(16*epsilon^3) * (9*x^3 + 27*epsilon*x^2 + 27*epsilon^2*x + 9*epsilon^3);
            elseif (x < epsilon/3)
                u = u_L/(16*epsilon^3) * (-18*x^3 + 18*epsilon^2*x + 8*epsilon^3);
            elseif (x < epsilon)
                u = u_L/(16*epsilon^3) * (9*x^3 - 27*epsilon*x^2 + 27*epsilon^2*x + 7*epsilon^3);
            else
                u = u_L;
            end
            
        case 8
            if (x < -epsilon)
                u = 0;
            elseif (x < -epsilon/2)
                u = u_L/(6*epsilon^4) * (4*x^4 + 16*epsilon*x^3 + 24*epsilon^2*x^2 + 16*epsilon^3*x + 4*epsilon^4);
            elseif (x < 0)
                u = u_L/(6*epsilon^4) * (-12*x^4 - 16*epsilon*x^3 + 8*epsilon^3*x + 3*epsilon^4);
            elseif (x < epsilon/2)
                u = u_L/(6*epsilon^4) * (12*x^4 - 16*epsilon*x^3 + 8*epsilon^3*x + 3*epsilon^4);
            elseif (x < epsilon)
                u = u_L/(6*epsilon^4) * (-4*x^4 + 16*epsilon*x^3 - 24*epsilon^2*x^2 + 16*epsilon^3*x + 2*epsilon^4);
            else
                u = u_L;
            end
            
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
