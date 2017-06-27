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
         ...
         'material_A'              , ...
         'material_G_c'            , ...
         'material_ell_0');
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'LM_array1', 'LM_array2'  , ...
         'BCU_array1', 'BCU_array2', ...
         'BCF_array1', 'BCF_array2', ...
         'index_u1', 'index_u2'    , ...
         'index_f1', 'index_f2'    , ...
         'u_L');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Initialize the displacement field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Load the patch assembly file for the displacement field
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
    %  Initialize the L2 projection matrix and the solution and RHS
    %  vectors
    %----------------------------------------------------------------------
    % Number of matrix entries that we compute per element
    numMatrixEntriesPerElement = numDOFsPerElement^2;
    
    % Row, column, value arrays for the L2 projection matrix K
    temp = zeros(numMatrixEntries(1), 1);
    
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
    
    clear temp;
    
    
    % Initialize the solution vector and prescribe the known coefficients
    u1 = zeros(numDOFs(1), 1);
    u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
    
    % Initialize the RHS vector and prescribe the known coefficients
    f = zeros(numDOFs(1), 1);
    f(BCF_array1(:, 1)) = BCF_array1(:, 2);
    
    
    % Initialization for element matrices K11, K12, K22
    elementMatrix = zeros(numDOFsPerElement);
    
    % Initialization for element vectors f1, f2
    elementVector = zeros(numDOFsPerElement, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(numDOFsPerElement, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : numMatrixEntriesPerElement)';
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over elements (e = e1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    lastIndex_for_K = 0;
    
    
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
        basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
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
        K_e = elementMatrix;
        f_e = elementVector;
        
        % Find the equation indices for the element
        index_equation1 = LM_array1(:, e);
        
        
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
            %  Evaluate the exact solution
            %--------------------------------------------------------------
            u_exact = function_u_exact(x(q_e), u_L, material_ell_0, order);
            
            
            %--------------------------------------------------------------
            %  Form the element matrix
            %--------------------------------------------------------------
            K_e = K_e + w(q_e) * Jacobian * (basis_physical_der0 * basis_physical_der0');
            
            
            %--------------------------------------------------------------
            %  Form the element RHS vector
            %--------------------------------------------------------------
            f_e = f_e + w(q_e) * Jacobian * (u_exact * basis_physical_der0);
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over quadrature points
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Add the element matrix
        index = lastIndex_for_K + indices_for_elementMatrix;
        
        rows_for_K(index)    = kron(vector_of_ones, index_equation1);
        columns_for_K(index) = kron(index_equation1, vector_of_ones);
        values_for_K(index)  = reshape(K_e, numMatrixEntriesPerElement, 1);
        
        lastIndex_for_K = lastIndex_for_K + numMatrixEntriesPerElement;
        
        
        % Add the element RHS vector
        f(index_equation1) = f(index_equation1) + f_e;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Solve for the displacement field
    %----------------------------------------------------------------------
    % Assemble the L2 projection matrix
    K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs(1), numDOFs(1));
    
    
    % Apply a row preconditioner
    precond = spdiags(1./max(K, [], 2), 0, numDOFs(1), numDOFs(1));
    K = precond * K;
    f = precond * f;
    
    % Solve for the unknown coefficients
    u1(index_u1) = K(index_u1, index_u1) \ (f(index_u1) - K(index_u1, index_f1) * u1(index_f1));
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Initialize the phase field
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Load the patch assembly file for the phase field
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 2), ...
         'nodes'             , ...
         'bezierExtractions1', ...
         'IEN_array');
    
    
    %----------------------------------------------------------------------
    %  Initialize the solution and RHS vectors
    %----------------------------------------------------------------------
    % Initialize the solution vector and prescribe the known coefficients
    u2 = zeros(numDOFs(2), 1);
    u2(BCU_array2(:, 1)) = BCU_array2(:, 2);
    
    % Initialize the RHS vector and prescribe the known coefficients
    f = zeros(numDOFs(2), 1);
    f(BCF_array2(:, 1)) = BCF_array2(:, 2);
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over elements (e = e1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    lastIndex_for_K = 0;
    
    
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
        basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
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
        K_e = elementMatrix;
        f_e = elementVector;
        
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
            basis_physical_der0 = basis_der0(:, q_e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the exact solution
            %--------------------------------------------------------------
            c_exact = function_c_exact(x(q_e), material_ell_0, order);
            
            
            %--------------------------------------------------------------
            %  Form the element matrix
            %--------------------------------------------------------------
            K_e = K_e + w(q_e) * Jacobian * (basis_physical_der0 * basis_physical_der0');
            
            
            %--------------------------------------------------------------
            %  Form the element RHS vector
            %--------------------------------------------------------------
            f_e = f_e + w(q_e) * Jacobian * (c_exact * basis_physical_der0);
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over quadrature points
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Add the element matrix
        index = lastIndex_for_K + indices_for_elementMatrix;
        
        rows_for_K(index)    = kron(vector_of_ones, index_equation2);
        columns_for_K(index) = kron(index_equation2, vector_of_ones);
        values_for_K(index)  = reshape(K_e, numMatrixEntriesPerElement, 1);
        
        lastIndex_for_K = lastIndex_for_K + numMatrixEntriesPerElement;
        
        
        % Add the element RHS vector
        f(index_equation2) = f(index_equation2) + f_e;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Solve for the displacement field
    %----------------------------------------------------------------------
    % Assemble the L2 projection matrix
    K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs(2), numDOFs(2));
    
    
    % Apply a row preconditioner
    precond = spdiags(1./max(K, [], 2), 0, numDOFs(2), numDOFs(2));
    K = precond * K;
    f = precond * f;
    
    % Solve for the unknown coefficients
    u2(index_u2) = K(index_u2, index_u2) \ (f(index_u2) - K(index_u2, index_f2) * u2(index_f2));
    
    
    
    
    %----------------------------------------------------------------------
    %  Save the fields
    %----------------------------------------------------------------------
    save(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', 'u2', '-v6');
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
                u = u_L/(2*epsilon^2) * ( x^2 + 2*epsilon*x + epsilon^2);
            elseif (x < epsilon)
                u = u_L/(2*epsilon^2) * (-x^2 + 2*epsilon*x + epsilon^2);
            else
                u = u_L;
            end
            
        case 6
            if (x < -epsilon)
                u = 0;
            elseif (x < -epsilon/3)
                u = u_L/(16*epsilon^3) * (  9*x^3 + 27*epsilon*x^2 + 27*epsilon^2*x + 9*epsilon^3);
            elseif (x < epsilon/3)
                u = u_L/(16*epsilon^3) * (-18*x^3                  + 18*epsilon^2*x + 8*epsilon^3);
            elseif (x < epsilon)
                u = u_L/(16*epsilon^3) * (  9*x^3 - 27*epsilon*x^2 + 27*epsilon^2*x + 7*epsilon^3);
            else
                u = u_L;
            end
            
        case 8
            if (x < -epsilon)
                u = 0;
            elseif (x < -epsilon/2)
                u = u_L/(6*epsilon^4) * (  4*x^4 + 16*epsilon*x^3 + 24*epsilon^2*x^2 + 16*epsilon^3*x + 4*epsilon^4);
            elseif (x < 0)
                u = u_L/(6*epsilon^4) * (-12*x^4 - 16*epsilon*x^3                    +  8*epsilon^3*x + 3*epsilon^4);
            elseif (x < epsilon/2)
                u = u_L/(6*epsilon^4) * ( 12*x^4 - 16*epsilon*x^3                    +  8*epsilon^3*x + 3*epsilon^4);
            elseif (x < epsilon)
                u = u_L/(6*epsilon^4) * ( -4*x^4 + 16*epsilon*x^3 - 24*epsilon^2*x^2 + 16*epsilon^3*x + 2*epsilon^4);
            else
                u = u_L;
            end
            
    end
end


function c = function_c_exact(x, material_ell_0, order)
    % Normalize the physical domain
    y = abs(x)/material_ell_0;
    
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
