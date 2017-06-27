%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 3D J2 plasticity for isotropic and
%  homogeneous materials. Linear isotropic and kinematic hardening can be
%  accomodated.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_3d.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      restartTime is the time step at which we start the simulation
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement fields (.mat files)
%--------------------------------------------------------------------------
function model_3d_parallel(path_to_assembly_directory, path_to_results_directory, restartTime)
    myCluster = parcluster('local');
    myCluster.NumWorkers = 4;
    saveProfile(myCluster);
    
    pool_handle = parpool(4);
    
    
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  3D J2 plasticity for isotropic, homogeneous materials.\n\n');
    
    
    % Load the global assembly file
    file_assembly_global = sprintf('%sfile_assembly_global.mat', path_to_assembly_directory);
    data1 = load(file_assembly_global, ...
         'numPatches'              , ...
         'numNodesBeforePatch'     , ...
         'numMatrixEntries'        , ...
         'numDOFs'                 , ...
         'numDOFsPerNode'          , ...
         'GN_array'                , ...
         ...
         'material_sigmaY'         , ...
         'material_H'              , ...
         'isotropicHardeningFactor', ...
         'numLoadSteps'            , ...
         'numStressPoints'         , ...
         'numTimeStepsBetweenSaves', ...
         'maxNewtonsMethod'        , ...
         'tolNewtonsMethod');
    
    
    % Load the patch assembly file
    file_assembly_patch = sprintf('%sfile_assembly_patch%d.mat', path_to_assembly_directory, 1);
    data2 = load(file_assembly_patch, ...
         'elementType'       , ...
         'material_E'        , ...
         'material_nu'       , ...
         'p1'                , ...
         'p2'                , ...
         'p3'                , ...
         'nodes'             , ...
         'numNodesPerElement', ...
         'numDOFsPerElement' , ...
         'bezierExtractions1', ...
         'bezierExtractions2', ...
         'bezierExtractions3', ...
         'numElements1'      , ...
         'numElements2'      , ...
         'numElements3'      , ...
         'elementSizes1'     , ...
         'elementSizes2'     , ...
         'elementSizes3'     , ...
         'IEN_array'         , ...
         'numQuadraturePoints');
    
    
    numPatches = data1.numPatches;
    numNodesBeforePatch = data1.numNodesBeforePatch;
    numMatrixEntries = data1.numMatrixEntries;
    numDOFs = data1.numDOFs;
    numDOFsPerNode = data1.numDOFsPerNode;
    GN_array = data1.GN_array;
    material_sigmaY = data1.material_sigmaY;
    material_H = data1.material_H;
    isotropicHardeningFactor = data1.isotropicHardeningFactor;
    numLoadSteps = data1.numLoadSteps;
    numStressPoints = data1.numStressPoints;
    numTimeStepsBetweenSaves = data1.numTimeStepsBetweenSaves;
    maxNewtonsMethod = data1.maxNewtonsMethod;
    tolNewtonsMethod = data1.tolNewtonsMethod;
    
    elementType = data2.elementType;
    material_E = data2.material_E;
    material_nu = data2.material_nu;
    p1 = data2.p1;
    p2 = data2.p2;
    p3 = data2.p3;
    nodes = data2.nodes;
    numNodesPerElement = data2.numNodesPerElement;
    numDOFsPerElement = data2.numDOFsPerElement;
    bezierExtractions1 = data2.bezierExtractions1;
    bezierExtractions2 = data2.bezierExtractions2;
    bezierExtractions3 = data2.bezierExtractions3;
    numElements1 = data2.numElements1;
    numElements2 = data2.numElements2;
    numElements3 = data2.numElements3;
    elementSizes1 = data2.elementSizes1;
    elementSizes2 = data2.elementSizes2;
    elementSizes3 = data2.elementSizes3;
    IEN_array = data2.IEN_array;
    numQuadraturePoints = data2.numQuadraturePoints;
    
    
    
    addAttachedFiles(pool_handle, {file_assembly_global, file_assembly_patch});
    
    
    % Read the element type
    switch (lower(elementType))
        case {'bspline', 'b-spline'}
            elementType = 0;
            
        case 'nurbs'
            elementType = 1;
            
        otherwise
            fprintf('\n');
            fprintf('  Error: Element type must be B-splines or NURBS. The problem will be left unsolved.\n\n');
            
            quit;
    end
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Pre-process
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    if (restartTime == 0)
        %------------------------------------------------------------------
        %  Initialize loading
        %------------------------------------------------------------------
        % Set the load step (m) and time step (n) to the beginning
        m = 1;
        n = 1;
        
        % Initialize the total displacements and rotations applied
        displacementSoFar = [0; 0; 0];
        rotationSoFar = 0;
        
        % Initialize the increment level
        incrementLevel = 0;
        
        % Initialize the number of successive convergences
        numSuccessiveConvergences = 0;
        
        % Flag that indicates whether the load step is finished
        flag_isLoaded = 0;
        
        
        %------------------------------------------------------------------
        %  Initialize the fields
        %  
        %  We use Voigt notation and view the 3 x 3 strain and stress
        %  tensors as 6 x 1 vectors. The entries for the (total) strain
        %  and the plastic strain vectors are as follows:
        %  
        %    strain = [    strain_11;
        %                  strain_22;
        %                  strain_33;
        %              2 * strain_23;
        %              2 * strain_13;
        %              2 * strain_12];    (similarly for plastic strain)
        %  
        %  The entries of the stress vector are as follows:
        %  
        %    stress = [stress_11;
        %              stress_22;
        %              stress_33;
        %              stress_23;
        %              stress_13;
        %              stress_12];
        %  
        %------------------------------------------------------------------
        temp = zeros(6, numStressPoints);
        
        % Assume zero displacements at t = 0
        u = zeros(numDOFs, 1);
        
        % Stress tensor at the current time step
%        stress = temp;
        
        % Backstress tensor at the current and previous time steps
%        backstress     = temp;
        backstress_old = temp;
        
        % Strain tensor at the current time step
%        strain = temp;
        
        % Plastic strain tensor at the current and previous time steps
%        strain_pl     = temp;
        strain_pl_old = temp;
        
        % Equivalent plastic strain at the current and previous time steps
%        strain_pl_eq     = zeros(1, numStressPoints);
        strain_pl_eq_old = zeros(1, numStressPoints);
        
        clear temp;
        
    else
        %------------------------------------------------------------------
        %  Initialize loading
        %------------------------------------------------------------------
        % Load the results from a previous time step
        load(sprintf('%sfile_results_time%06.0f', path_to_results_directory, restartTime), ...
             'u'                        , ...
             'stress'                   , ...
             'backstress'               , ...
             'strain'                   , ...
             'strain_pl'                , ...
             'strain_pl_eq'             , ...
             ...
             'loadStep'                 , ...
             'displacementSoFar'        , ...
             'rotationSoFar'            , ...
             'incrementLevel'           , ...
             'numSuccessiveConvergences', ...
             'flag_isLoaded'            , ...
             '-v7.3');
        
        % Set the load step (m) and time step (n)
        m = loadStep;
        n = restartTime + 1;
        
        
        %------------------------------------------------------------------
        %  Initialize the fields
        %------------------------------------------------------------------
        % Read the backstress tensor, plastic strain tensor, and equivalent
        % plastic strain at the previous time step
        backstress_old = backstress;
        strain_pl_old = strain_pl;
        strain_pl_eq_old = strain_pl_eq;
        
        clear backstress strain_pl strain_pl_eq;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants for quadrature
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    constant_p3p1 = p3 + 1;
    
    % Set the quadrature rule
    [z1, z2, z3, w] = set_3d_gauss_quadrature_for_bernstein(numQuadraturePoints);
    numQuadraturePointsPerElement = prod(numQuadraturePoints);
    
    
    % Evaluate the Bernstein polynomials for direction 1 at
    % the quadrature points
    Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints(1));
    Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints(1));
    
    for j = 1 : numQuadraturePoints(1)
        temp = eval_1d_bernstein_der(z1(j), p1);
        
        Bernstein1_der0(:, j) = temp(:, 1);
        Bernstein1_der1(:, j) = temp(:, 2);
    end
    
    % Evaluate the Bernstein polynomials for direction 2 at
    % the quadrature points
    Bernstein2_der0 = zeros(constant_p2p1, numQuadraturePoints(2));
    Bernstein2_der1 = zeros(constant_p2p1, numQuadraturePoints(2));
    
    for j = 1 : numQuadraturePoints(2)
        temp = eval_1d_bernstein_der(z2(j), p2);
        
        Bernstein2_der0(:, j) = temp(:, 1);
        Bernstein2_der1(:, j) = temp(:, 2);
    end
    
    % Evaluate the Bernstein polynomials for direction 3 at
    % the quadrature points
    Bernstein3_der0 = zeros(constant_p3p1, numQuadraturePoints(3));
    Bernstein3_der1 = zeros(constant_p3p1, numQuadraturePoints(3));
    
    for j = 1 : numQuadraturePoints(3)
        temp = eval_1d_bernstein_der(z3(j), p3);
        
        Bernstein3_der0(:, j) = temp(:, 1);
        Bernstein3_der1(:, j) = temp(:, 2);
    end
    
    clear z1 z2 z3 temp;
    
    
    %----------------------------------------------------------------------
    %  Build the elasticity matrix
    %----------------------------------------------------------------------
    % Bulk modulus
    material_K = material_E / (3 * (1 - 2 * material_nu));
    
    % Shear modulus
    material_G = material_E / (2 * (1 + material_nu));
    
    % D11, D22, D33 entries
    D11 = material_K + 4/3 * material_G;
    
    % D12, D13, D23, D21, D31, D32 entries
    D12 = material_K - 2/3 * material_G;
    
    % D44, D55, D66 entries
    D33 = material_G;
    
    % Set the elasticity matrix
    D_el = [D11, D12, D12,   0,   0,   0; ...
            D12, D11, D12,   0,   0,   0; ...
            D12, D12, D11,   0,   0,   0; ...
              0,   0,   0, D33,   0,   0; ...
              0,   0,   0,   0, D33,   0; ...
              0,   0,   0,   0,   0, D33];
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Pre-process
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over load steps (m)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for plasticity
    constant_sq2d3 = sqrt(2/3);
    constant_2G = 2 * material_G;
    
    
    for loadStep = m : numLoadSteps
        % Find the positions of the nodes at the beginning of the load step
        % (TODO)
        nodesAtBeginning = nodes(:, [1; 2; 3])';
        
        
        % Load the load step file
        load(sprintf('%sfile_loadstep%d', path_to_assembly_directory, loadStep), ...
             'numBoundaries', ...
             'numTimeSteps');
        
        % Initialize the increment
        increment = (1 / numTimeSteps) * 2^incrementLevel;
        
        % Initialize the ID array
        ID_array = [];
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over time steps (n)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Continue to loop until we reach the total displacement or the
        % total rotation that is specified by the load step
        while (flag_isLoaded ~= 1)
            fprintf('\n');
            fprintf('- Load step = %d, time step = %d\n', loadStep, n);
            fprintf('  Increment level = %d\n', incrementLevel);
            
            
            %--------------------------------------------------------------
            %  Set the displacement and force BCs
            %--------------------------------------------------------------
            BCs_displacement = [];
            BCs_force = [];
            
            for boundary = 1 : numBoundaries
                % Load the BC file
                load(sprintf('%sfile_loadstep%d_boundary%d', path_to_assembly_directory, loadStep, boundary), ...
                     'loadType', ...
                     'nodesOnBoundary', ...
                     'dofIndex', ...
                     'displacementTotal', ...
                     'rotationTotal', ...
                     'rotationAxis');
                
                numNodesOnBoundary = size(nodesOnBoundary, 1);
                numDOFIndex = size(dofIndex, 1);
                
                % Displacement only
                if (loadType == 1)
                    % Evaluate the incremental displacements for the nodes
                    % on the boundary
                    bcValue = increment * displacementTotal;
                    
                    BCs_displacement = [BCs_displacement; ...
                                        [kron(nodesOnBoundary, ones(numDOFIndex, 1)), repmat(dofIndex, numNodesOnBoundary, 1), repmat(bcValue, numNodesOnBoundary, 1)]];
                    
                % Rotation only
                elseif (loadType == 2)
                    % Evaluate the rotation matrix that maps the nodes on 
                    % the boundary at the beginning of the load step to
                    % those at the current time step
                    theta0 = rotationSoFar;
                    
                    constant_costheta = cos(theta0);
                    constant_sintheta = sin(theta0);
                    constant_1mcostheta = 1 - constant_costheta;
                    
                    rotationMatrix0 = constant_1mcostheta * (rotationAxis * rotationAxis') + ...
                                      [                   constant_costheta, -rotationAxis(3) * constant_sintheta,  rotationAxis(2) * constant_sintheta; ...
                                        rotationAxis(3) * constant_sintheta,                    constant_costheta, -rotationAxis(1) * constant_sintheta; ...
                                       -rotationAxis(2) * constant_sintheta,  rotationAxis(1) * constant_sintheta,                    constant_costheta];
                    
                    % Evaluate the rotation matrix that maps the nodes on 
                    % the boundary at the beginning of the load step to
                    % those at the next time step
                    theta1 = rotationSoFar + increment * rotationTotal;
                    
                    constant_costheta = cos(theta1);
                    constant_sintheta = sin(theta1);
                    constant_1mcostheta = 1 - constant_costheta;
                    
                    rotationMatrix1 = constant_1mcostheta * (rotationAxis * rotationAxis') + ...
                                      [                   constant_costheta, -rotationAxis(3) * constant_sintheta,  rotationAxis(2) * constant_sintheta; ...
                                        rotationAxis(3) * constant_sintheta,                    constant_costheta, -rotationAxis(1) * constant_sintheta; ...
                                       -rotationAxis(2) * constant_sintheta,  rotationAxis(1) * constant_sintheta,                    constant_costheta];
                    
                    
                    % Evaluate the incremental displacements for the nodes
                    % on the boundary
                    bcValue = (rotationMatrix1 - rotationMatrix0) * nodesAtBeginning(dofIndex, nodesOnBoundary);
                    
                    BCs_displacement = [BCs_displacement; ...
                                        [kron(nodesOnBoundary, ones(numDOFIndex, 1)), repmat(dofIndex, numNodesOnBoundary, 1), reshape(bcValue, numNodesOnBoundary * numDOFIndex, 1)]];
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Build the ID, LM, and BC arrays
            %--------------------------------------------------------------
            % Unless the nodes on the boundaries change within a load step,
            % we only need to build the ID and LM arrays once
            if (isempty(ID_array))
                ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array);
                LM_array = build_lm_array(IEN_array, ID_array, GN_array, 0);
            end
            
            % We build the BC arrays at every time step since the BC value
            % differ from one time step to another. The equation indices
            % will stay the same, however.
            [BCU_array, BCF_array] = build_bc_array(BCs_displacement, BCs_force, ID_array);
            
            % Number of unknown degrees of freedom
            numUnknownDOFs = numDOFs - size(BCU_array, 1);
            
            % Indices for unknown displacements
            index_u = (1 : numUnknownDOFs)';
            index_f = ((numUnknownDOFs + 1) : numDOFs)';
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over Newton's method (k)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            %  Initialize the row, column, and value arrays for the consistent
            %  tangent matrix
            %--------------------------------------------------------------
            % By definition, the initial Newton's guess is the solution
            % vector from the previous time step
            u_new = u;
            
            % Initialize the Newton's increment
            u_increment = zeros(numDOFs, 1);
            
            % External force vector due to increments in body force and
            % tractions
            f_external = zeros(numDOFs, 1);
            f_external(BCF_array(:, 1)) = BCF_array(:, 2);
            
            % Flag that indicates whether we are close to convergence
            flag_isCloseToConvergence = 0;
            
            % Norm of the residual vector
            residual_norm = 0;
            
            
            for k = 0 : maxNewtonsMethod
                rows_for_K    = [];
                columns_for_K = [];
                values_for_K  = [];
                
                rows_for_f    = [];
                columns_for_f = [];
                values_for_f  = [];
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over elements (e)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Counters for the element and the quadrature point
                numElements = numElements1 * numElements2 * numElements3;
                
                backstress = [];
                stress = [];
                strain = [];
                strain_pl = [];
                strain_pl_eq = [];
                flag_isPlastic = [];
                
                parfor e = 1 : numElements
                    e1 = mod(e - 1, numElements1) + 1;
                    e2 = mod(ceil(e / numElements1) - 1 , numElements2) + 1;
                    e3 = mod(ceil(e / (numElements1 * numElements2)) - 1, numElements3) + 1;
                    
                    
                    % Find the Bezier extraction matrix
                    bezierExtractions1_e = bezierExtractions1(:, :, e1);
                    
                    % Evaluate the map derivative dt/dxi (constant)
                    dt_dxi = 1 / elementSizes1(e1);
                    
                    % Evaluate the univariate B-splines
                    Bspline1_der0 = bezierExtractions1_e * Bernstein1_der0;
                    Bspline1_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
                    
                    
                    % Find the Bezier extraction matrix
                    bezierExtractions2_e = bezierExtractions2(:, :, e2);
                    
                    % Evaluate the map derivative dt/deta (constant)
                    dt_deta = 1 / elementSizes2(e2);
                    
                    % Evaluate the univariate B-splines
                    Bspline2_der0 = bezierExtractions2_e * Bernstein2_der0;
                    Bspline2_der1 = (dt_deta * bezierExtractions2_e) * Bernstein2_der1;
                    
                    
                    % Find the Bezier extraction matrix
                    bezierExtractions3_e = bezierExtractions3(:, :, e3);
                    
                    % Evaluate the map derivative dt/dzeta (constant)
                    dt_dzeta = 1 / elementSizes3(e3);
                    
                    % Evaluate the univariate B-splines
                    Bspline3_der0 = bezierExtractions3_e * Bernstein3_der0;
                    Bspline3_der1 = (dt_dzeta * bezierExtractions3_e) * Bernstein3_der1;
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the basis functions in the parametric domain
                    %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                    %------------------------------------------------------
                    % For B-splines
                    if (elementType == 0)
                        % Find the positions of the nodes
                        nodes_e = nodes(IEN_array(:, e), :)';
                        
                        % Evaluate the B-splines
                        basis_der100 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der1));
                        basis_der010 = kron(Bspline3_der0, kron(Bspline2_der1, Bspline1_der0));
                        basis_der001 = kron(Bspline3_der1, kron(Bspline2_der0, Bspline1_der0));
                        
                    % For NURBS
                    elseif (elementType == 1)
                        % Find the positions and weights of the nodes
                        nodes_e = nodes(IEN_array(:, e), [1; 2; 3])';
                        w_e = nodes(IEN_array(:, e), 4);
                        
                        % Evaluate the B-splines
                        Bspline_der000 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der0));
                        Bspline_der100 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der1));
                        Bspline_der010 = kron(Bspline3_der0, kron(Bspline2_der1, Bspline1_der0));
                        Bspline_der001 = kron(Bspline3_der1, kron(Bspline2_der0, Bspline1_der0));
                        
                        % Evaluate the NURBS
                        [~, basis_der100, basis_der010, basis_der001] = eval_3d_nurbs_der1(w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001);
                        
                    end
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the map derivatives (x = [x1; x2; x3])
                    %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                    %------------------------------------------------------
                    dx_dxi = nodes_e * basis_der100;
                    dx_deta = nodes_e * basis_der010;
                    dx_dzeta = nodes_e * basis_der001;
                    
                    
                    
                    %------------------------------------------------------
                    % -----------------------------------------------------
                    %   Begin: Loop over quadrature points (q, q_e)
                    % -----------------------------------------------------
                    %------------------------------------------------------
                    % Initialize the element consistent tangent matrix and
                    % the element internal force vector
                    K_e = zeros(numDOFsPerElement);
                    f_e = zeros(numDOFsPerElement, 1);
                    
                    % Get the current Newton's guess for the solution vector
                    index_equation = LM_array(:, e);
                    u_e = u_new(index_equation);
                    
                    
                    for q_e = 1 : numQuadraturePointsPerElement
                        q = numQuadraturePointsPerElement * (e - 1) + q_e;
                        
                        %--------------------------------------------------
                        %  Form the Jacobian matrix
                        %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                        %--------------------------------------------------
                        JacobianMatrix = [dx_dxi(:, q_e), dx_deta(:, q_e), dx_dzeta(:, q_e)]';
                        
                        % Evaluate the Jacobian
                        Jacobian = det(JacobianMatrix) / (dt_dxi * dt_deta * dt_dzeta);
                        
                        if (Jacobian <= 0)
                            fprintf('\n');
                            fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                            
                            quit;
                        end
                        
                        
                        %--------------------------------------------------
                        %  Evaluate the basis functions in the physical domain
                        %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                        %--------------------------------------------------
                        basis_der1 = [basis_der100(:, q_e), basis_der010(:, q_e), basis_der001(:, q_e)]';
                        
                        basis_physical_der1 = JacobianMatrix \ basis_der1;
                           
                        
                        %--------------------------------------------------
                        %  Evaluate the element B matrix, which gives the
                        %  element strains when multiplied by the element
                        %  displacement vector
                        %--------------------------------------------------
                        % Initialize the B matrix
                        B = zeros(6, numDOFsPerElement);
                        
                        % Column indices for the B matrix
                        temp1 = 1;
                        temp2 = 2;
                        temp3 = 3;
                        
                        for a = 1 : numNodesPerElement
                            B(1, temp1) = basis_physical_der1(1, a);
                            B(5, temp1) = basis_physical_der1(3, a);
                            B(6, temp1) = basis_physical_der1(2, a);
                            
                            B(2, temp2) = basis_physical_der1(2, a);
                            B(4, temp2) = basis_physical_der1(3, a);
                            B(6, temp2) = basis_physical_der1(1, a);
                            
                            B(3, temp3) = basis_physical_der1(3, a);
                            B(4, temp3) = basis_physical_der1(2, a);
                            B(5, temp3) = basis_physical_der1(1, a);
                            
                            temp1 = temp1 + 3;
                            temp2 = temp2 + 3;
                            temp3 = temp3 + 3;
                        end
                        
                        
                        %--------------------------------------------------
                        %  Evaluate the trial strains and stresses
                        %--------------------------------------------------
                        % Evaluate the strain tensor
                        strain_q = B * u_e;
                        
                        % Evaluate the elastic strain tensor and its trace
                        strain_el = strain_q - strain_pl_old(:, q);
                        strain_el_tr = strain_el(1) + strain_el(2) + strain_el(3);
                        
                        % Evaluate the stress tensor and its trace
                        stress_q = D_el * strain_el;
                        stress_tr = stress_q(1) + stress_q(2) + stress_q(3);
                        
                        % Evaluate the deviatoric stress tensor
                        stress_dev = stress_q;
                        stress_dev(1) = stress_dev(1) - stress_tr / 3;
                        stress_dev(2) = stress_dev(2) - stress_tr / 3;
                        stress_dev(3) = stress_dev(3) - stress_tr / 3;
                        
                        % Evaluate the deviatoric stress tensor,
                        % relative to the backstress tensor
                        stress_rel = stress_dev - backstress_old(:, q);
                        
                        
                        %--------------------------------------------------
                        %  Check the J2 yield condition
                        %--------------------------------------------------
                        % Evaluate the norm of the deviatoric stress tensor
                        stress_dev_norm = sqrt(stress_rel(1)^2 + stress_rel(2)^2 + stress_rel(3)^2 + 2 * (stress_rel(4)^2 + stress_rel(5)^2 + stress_rel(6)^2));
                        
                        % Evaluate the hardening function
                        hardening_function = material_sigmaY + isotropicHardeningFactor * material_H * strain_pl_eq_old(q);
                        
                        % Evaluate the yield function
                        yield_function = stress_dev_norm - constant_sq2d3 * hardening_function;
                        
                        
                        %--------------------------------------------------
                        %  Elastic state
                        %--------------------------------------------------
                        if (yield_function <= 0)
                            % Set the flag to 0
                            flag_isPlastic_q = 0;
                            
                            
                            %----------------------------------------------
                            %  Accept the trial state
                            %----------------------------------------------
                            % Do not update the backstress tensor
                            backstress_q = backstress_old(:, q);
                            
                            % Do not update the plastic strain tensor
                            strain_pl_q = strain_pl_old(:, q);
                            
                            % Do not update the equivalent plastic strain
                            strain_pl_eq_q = strain_pl_eq_old(q);
                            
                            
                            %----------------------------------------------
                            %  Evaluate the consistent tangent matrix and
                            %  the internal force vector
                            %----------------------------------------------
                            K_e = K_e + B' * (w(q_e) * Jacobian * D_el) * B;
                            f_e = f_e + B' * (w(q_e) * Jacobian * stress_q);
                            
                            
                        %----------------------------------------------
                        %  Plastic state
                        %----------------------------------------------
                        else
                            % Set the flag to 1
                            flag_isPlastic_q = 1;
                            
                            
                            %----------------------------------------------
                            %  Evaluate the consistency parameter
                            %----------------------------------------------
                            % For nonlinear isotropic and kinematic hardening, 
                            % we would need to use Newton's method instead
                            Delta_gamma = yield_function / (constant_2G + 2/3 * material_H);
                            
                            % Some useful constant for radial return
                            constant1 = Delta_gamma * constant_2G;
                            
                            
                            %----------------------------------------------
                            %  Radially return to the yield surface
                            %----------------------------------------------
                            % Evaluate the normal tensor, which is the unit 
                            % deviatoric stress tensor
                            normal = stress_rel / stress_dev_norm;
                            
                            % Update the backstress tensor
                            backstress_q = backstress_old(:, q) + (constant_sq2d3 * (1 - isotropicHardeningFactor) * material_H) * normal;
                            
                            % Update the stress tensor
                            stress_q = stress_dev - constant1 * normal;
                            stress_q(1) = stress_q(1) + material_K * strain_el_tr;
                            stress_q(2) = stress_q(2) + material_K * strain_el_tr;
                            stress_q(3) = stress_q(3) + material_K * strain_el_tr;
                            
                            % Update the plastic strain tensor
                            strain_pl_q = strain_pl_old(:, q) + Delta_gamma * [    normal(1); ...
                                                                                   normal(2); ...
                                                                                   normal(3); ...
                                                                               2 * normal(4); ...
                                                                               2 * normal(5); ...
                                                                               2 * normal(6)];
                            
                            % Update the equivalent plastic strain
                            strain_pl_eq_q = strain_pl_eq_old(q) + Delta_gamma * constant_sq2d3;
                            
                            
                            %----------------------------------------------
                            %  Evaluate the consistent tangent matrix and
                            %  the internal force vector
                            %----------------------------------------------
                            % Some useful constants for the consistent tangent matrix
                            % (note, we multiplied theta and theta_bar by 2G)
                            constant2 = constant_2G * (1 - constant1 / stress_dev_norm);
                            constant3 = constant_2G * (1 - constant2 / constant_2G - material_G / (material_G + material_H / 3));
                            
                            % D11, D22, D33 entries
                            D11 = material_K + 2/3 * constant2;
                            
                            % D12, D13, D23, D21, D31, D32 entries
                            D12 = material_K - 1/3 * constant2;
                            
                            % D44, D55, D66 entries
                            D33 = constant2 / 2;
                            
                            % Set the elastoplastic matrix
                            D_ep = [D11, D12, D12,   0,   0,   0; ...
                                    D12, D11, D12,   0,   0,   0; ...
                                    D12, D12, D11,   0,   0,   0; ...
                                      0,   0,   0, D33,   0,   0; ...
                                      0,   0,   0,   0, D33,   0; ...
                                      0,   0,   0,   0,   0, D33] ...
                                 + constant3 * (normal * normal');
                            
                            K_e = K_e + B' * (w(q_e) * Jacobian * D_ep) * B;
                            f_e = f_e + B' * (w(q_e) * Jacobian * stress_q);
                            
                        end
                        
                        
                        stress = [stress, stress_q];
                        backstress = [backstress, backstress_q];
                        strain = [strain, strain_q];
                        strain_pl = [strain_pl, strain_pl_q];
                        strain_pl_eq = [strain_pl_eq, strain_pl_eq_q];
                        flag_isPlastic = [flag_isPlastic, flag_isPlastic_q];
                    end
                    
                    
                    %------------------------------------------------------
                    % -----------------------------------------------------
                    %   End: Loop over quadrature points
                    % -----------------------------------------------------
                    %------------------------------------------------------
                    
                    
                    %------------------------------------------------------
                    %  Global assembly
                    %------------------------------------------------------
                    temp = ones(numDOFsPerElement, 1);
                    
                    % Add the element matrix
                    for j = 1 : numDOFsPerElement
                        rows_for_K    = [rows_for_K; index_equation];
                        columns_for_K = [columns_for_K; index_equation(j) * temp];
                        values_for_K  = [values_for_K; K_e(:, j)];
                    end
                    
                    % Add the element RHS vector
                    rows_for_f    = [rows_for_f; index_equation];
                    columns_for_f = [columns_for_f; temp];
                    values_for_f  = [values_for_f; f_e];
                    
                end
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Loop over elements
                % ---------------------------------------------------------
                %----------------------------------------------------------
                
                
                
                %----------------------------------------------------------
                %  Check for convergence
                %  
                %  We accept that the Newton's guess corresponds to an
                %  equilibrium state via a two-step check.
                %----------------------------------------------------------
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
                        
                        % Set the flag to 1
                        flag_hasConverged = 1;
                        
                        % Record the number of successive convergences
                        if (incrementLevel < 0)
                            numSuccessiveConvergences = numSuccessiveConvergences + 1;
                        end
                        
                        break;
                        
                    elseif (k == maxNewtonsMethod)
                        fprintf('\n');
                        fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                        
                        % Set the flag to 0
                        flag_hasConverged = 0;
                        
                        % Reset the number of successive convergences
                        numSuccessiveConvergences = 0;
                        
                        break;
                        
                    end
                end
                
                
                %----------------------------------------------------------
                %  Solve for the Newton's increment
                %----------------------------------------------------------
                % Assemble the consistent tangent matrix
                K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs, numDOFs);                
                f_internal = sparse(rows_for_f, columns_for_f, values_for_f, numDOFs, 1);
                
                % Save the norm of the residual vector at the previous iteration
                residual_norm0 = residual_norm;
                
                % Form the residual vector (note the negative sign)
                residual = f_external - f_internal;
                
                % Evaluate the norm of the residual vector at the current iteration
                residual_norm = norm(residual(index_u), 2);
                
                
                % Apply a row preconditioner
                precond = spdiags(1./max(K, [], 2), 0, numDOFs, numDOFs);
                K = precond * K;
                residual = precond * residual;
                
                
                % At the first iteration, we tell the Newton's increment to
                % satisfy any incremental displacements on the boundary
                if (k == 0)
                    u_increment(BCU_array(:, 1)) = BCU_array(:, 2);
                    u_increment(index_u) = K(index_u, index_u) \ (residual(index_u) - K(index_u, index_f) * u_increment(index_f));
                    
                % At subsequent iterations, the Newton's guess already satisfies
                % the displacement BCs. Hence, we tell the Newton's increment to
                % stop imposing the incremental displacements on the bodundary
                else
                    if (k == 1)
                        u_increment(BCU_array(:, 1)) = 0;
                    end
                    
                    u_increment(index_u) = K(index_u, index_u) \ residual(index_u);
                    
                end
                
                % Update the Newton's guess
                u_new = u_new + u_increment;
                
                
                % Evaluate the norm of the Newton's increment and check if
                % we are close to convergence
                increment_norm = norm(u_increment(index_u), 2);
                
                if (increment_norm < tolNewtonsMethod)
                    flag_isCloseToConvergence = 1;
                end
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over Newton's method
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Decide what to do next
            %--------------------------------------------------------------
            if (flag_hasConverged == 1)
                if (loadType == 1)
                    % Update the total displacement applied
                    displacementSoFar = displacementSoFar + increment * displacementTotal;
                    
                    % Check if we are close to finishing the load step
                    if (displacementSoFar(1) >= displacementTotal(1) && displacementSoFar(2) >= displacementTotal(2) && displacementSoFar(3) >= displacementTotal(3))
                        flag_isLoaded = 1;
                    end
                    
                elseif (loadType == 2)
                    % Update the total rotation applied
                    rotationSoFar = rotationSoFar + increment * rotationTotal;
                    
                    % Check if we are close to finishing the load step
                    if (rotationSoFar >= rotationTotal)
                        flag_isLoaded = 1;
                    end
                    
                end
                
                
                % Save the results
                if (mod(n, numTimeStepsBetweenSaves) == 0 || flag_isLoaded == 1)
                    save(sprintf('%sfile_results_time%06.0f', path_to_results_directory, n), ...
                         'ID_array', ...
                         'u', 'stress', 'backstress', 'strain', 'strain_pl', 'strain_pl_eq', 'flag_isPlastic', ...
                         'loadStep', 'displacementSoFar', 'rotationSoFar', 'incrementLevel', 'numSuccessiveConvergences', 'flag_isLoaded', ...
                         '-v7.3');
                end
                
                
                % Update the time step
                n = n + 1;
                
                % Update the fields
                u = u_new;
                backstress_old = backstress;
                strain_pl_old = strain_pl;
                strain_pl_eq_old = strain_pl_eq;
                
                % Check whether to relax the increment
                if (incrementLevel < 0 && numSuccessiveConvergences == 4)
                    increment = 2 * increment;
                    incrementLevel = incrementLevel + 1;
                end
                
            else
                % Lower the increment
                increment = 0.5 * increment;
                incrementLevel = incrementLevel - 1;
                
                % If the increment level is deemed too low, then terminate
                % and consider the problem unsolved
                if (incrementLevel < -7)
                    fprintf('\n');
                    fprintf('  Error: Equilibrium could not be reached with the lowest increment. The problem will be left unsolved.\n\n');
                    
                    quit;
                end
                
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over time step
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        % Reset the values for the next load step
        displacementSoFar = [0; 0; 0];
        rotationSoFar = 0;
        incrementLevel = 0;
        numSuccessiveConvergences = 0;
        flag_isLoaded = 0;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over load step
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
