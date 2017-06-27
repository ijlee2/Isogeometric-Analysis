%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the momentum equation for a 1D cracked bar in
%  tension.
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
%      order is the order of the phase field theory (2, 4, 6, 8)
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function model_1d(order, path_to_assembly_directory, path_to_results_directory)
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
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'LM_array1'      , 'LM_array2'      , ...
         'BCU_array1'     , 'BCU_array2'     , ...
         'BCF_array1'     , 'BCF_array2'     , ...
         'numUnknownDOFs1', 'numUnknownDOFs2', ...
         'index_u1'       , 'index_u2'       , ...
         'index_f1'       , 'index_f2');
    
    
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
    %   Begin: Find the displacement field (u1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for phase field theory
    constant_stiffness = material_E * material_A;
    
    
    % Some useful constant for global assembly
    constant_p1p1_sq = constant_p1p1^2;
    
    
    % Initialization for solution and RHS vectors
    globalVector  = zeros(numDOFs, 1);
    
    % Initialization for element matrices
    elementMatrix = zeros(constant_p1p1);
    
    % Initialization for element vectors
    elementVector = zeros(constant_p1p1, 1);
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(constant_p1p1, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : constant_p1p1_sq)';
    
    
    % Number of matrix entries that we compute for the momentum equation
    numMatrixEntriesPerElement = constant_p1p1_sq;
    
    % Initialize row, column, value arrays for the stiffness matrix K
    temp = zeros(numMatrixEntriesPerElement * numElements1, 1);
    
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
    
    clear temp;
    
    
    
    % Initialize the RHS vector and prescribe the known forces
    u1 = globalVector;
    f1 = globalVector;
    f1(BCF_array1(:, 1)) = BCF_array1(:, 2);
    
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
        basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
        basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
        
        
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
        % Initialize the element stiffness matrix and the element internal
        % force vector
        K11_e = elementMatrix;
%       f1_e  = elementVector;
        
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
            JacobianMatrix_inv = 1 / JacobianMatrix;
            
            basis_physical_der0 =                      basis_der0(:, q_e);
            basis_physical_der1 = JacobianMatrix_inv * basis_der1(:, q_e);
            
            
            %--------------------------------------------------------------
            %  Find the point on the physical domain
            %--------------------------------------------------------------
            x = nodes_e * basis_physical_der0;
            
            
            %--------------------------------------------------------------
            %  Evaluate the phase field
            %--------------------------------------------------------------
            y = abs(x) / material_ell_0;
            
            switch order
                case 2
                    c_der0 = 1 - exp(-1/2*y);
                    
                case 4
                    c_der0 = 1 - exp(-y) * (1 + y);
                    
                case 6
                    c_der0 = 1 - exp(-3/2*y) * (1 + 3/2*y + 9/8*y^2);
                    
                case 8
                    c_der0 = 1 - exp(-2*y) * (1 + 2*y + 2*y^2 + 4/3*y^3);
                    
            end
            
            
            %--------------------------------------------------------------
            %  Form the element matrix K11
            %--------------------------------------------------------------
            % Coefficients for the (regular) phase field theory.
            % The subscripts denote the derivatives of N_A and N_B.
            C11 = constant_stiffness * c_der0^2;
            
            K11_e = K11_e + (w(q_e) * Jacobian * C11) * (basis_physical_der1 * basis_physical_der1');
            
            
            %--------------------------------------------------------------
            %  Form the element RHS vector f1
            %--------------------------------------------------------------
%           f1_e = f1_e + (w(q_e) * Jacobian * 0) * basis_physical_der0;
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over quadrature points
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        %  Global assembly
        %------------------------------------------------------------------
        % Add the element matrix K11
        index = lastIndex_for_K + indices_for_elementMatrix;
        
        rows_for_K(index)    = kron(vector_of_ones, index_equation1);
        columns_for_K(index) = kron(index_equation1, vector_of_ones);
        values_for_K(index)  = reshape(K11_e, constant_p1p1_sq, 1);
        
        lastIndex_for_K = lastIndex_for_K + constant_p1p1_sq;
        
        
        % Add the element RHS vector f1
%       f1(index_equation1) = f1(index_equation1) + f1_e;
        
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over elements
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Solve for the displacement field
    %----------------------------------------------------------------------
    % Assemble the stiffness matrix
    K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs, numDOFs);
    
    
    % Apply a row preconditioner
    precond = spdiags(1./max(K, [], 2), 0, numDOFs, numDOFs);
    K  = precond * K;
    f1 = precond * f1;
    
    
    % Prescribe the known coefficients and solve for the unknown
    u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
    u1(index_u1) = K(index_u1, index_u1) \ (f1(index_u1) - K(index_u1, index_f1) * u1(index_f1));
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Save the results
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    save(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', '-v6');
end
