%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine creates the assembly file for the problem of 2D plate with
%  a circular hole in tension. We use quadratic NURBS as basis functions.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_generate_assemblyFile.m output
%  
%  Note that,
%  
%      numRefinements is the number of times we h-refine the knot vectors
%  
%  
%  Output:
%  
%  1. Assembly files (.mat files)
%--------------------------------------------------------------------------
function generate_assemblyFile_tractions(numRefinements)
    % Set the path to the assembly files directory
    directory_assembly = strcat('../assembly_files/');
    
    % Create the directory if it does not exist
    if ~exist(directory_assembly, 'dir')
        mkdir(directory_assembly);
    end
    
    % File path
    file_assembly = strcat(directory_assembly, 'file_assembly_numRefinements', sprintf('%d', numRefinements), '.mat');
    
    
    %----------------------------------------------------------------------
    %  Parameters for the material
    %----------------------------------------------------------------------
    % Material parameters
    material_E = 210e9;               % Young's modulus (Pa)
    material_nu = 0.3;                % Poisson's ratio (dimensionless)
    material_thickness = 0.001;       % Plate thickness (m)
    
    % Specimen size parameters
    plate_length = 0.1;               % The plate has length and width of 2L (m)
    hole_radius = 0.01;               % Radius of the circular hole (m)
    
    
    %----------------------------------------------------------------------
    %  Parameters for FEM
    %----------------------------------------------------------------------
    % Number of Gauss quadrature points
    numQuadraturePoints = [4; 4];
    
    
    %----------------------------------------------------------------------
    %  Create an array of NURBS structs (two patches)
    %  
    %  description
    %      -- the description of the NURBS patch used for debugging 
    %      -- string
    %  
    %  knots1, knots2, knots3
    %      -- the knot vectors in directions 1, 2, and 3
    %      -- column vectors
    %  
    %  nodes
    %      -- the control points in the physical space
    %      -- matrix with d columns, where d is the dimension
    %  
    %  p1, p2, p3
    %      -- the degree of the NURBS in directions 1, 2, and 3
    %      -- positive integers
    %----------------------------------------------------------------------
    %{
    numPatches = 1;
    
    patches(numPatches) = ...
        struct('type'               , 'nurbs', ...
               'description'        , '', ...
               'knots1'             , [], ...
               'knots2'             , [], ...
               'knots3'             , [], ...
               'nodes'              , [], ...
               'p1'                 , 0, ...
               'p2'                 , 0, ...
               'p3'                 , 0, ...
               'bezierExtractions1' , [], ...
               'bezierExtractions2' , [], ...
               'bezierExtractions3' , [], ...
               'numElements1'       , [], ...
               'numElements2'       , [], ...
               'numElements3'       , [], ...
               'elementSizes1'      , [], ...
               'elementSizes2'      , [], ...
               'elementSizes3'      , [], ...
               'sharedNodesIndex'   , [], ...
               'numDOFsPerNode'     , 0, ...
               'numDOFsPerElement'  , 0, ...
               'numDOFs'            , 0, ...
               'numUnknownDOFs'     , 0, ...
               'IEN_array'          , [], ...
               'ID_array'           , [], ...
               'LM_array'           , [], ...
               'BCU_array'          , [], ...
               'BCF_array'          , []);
    %}
    
    
    %----------------------------------------------------------------------
    %  Build patch 1
    %----------------------------------------------------------------------
    description = strcat('2D plate with a circular hole, refinement level', sprintf(' %d', numRefinements));
    
    numDOFsPerNode = 2;
    
    
    %----------------------------------------------------------------------
    %  Set the initial knots array
    %----------------------------------------------------------------------
    knots1 = [0; 0; 0; 1; 2; 2; 2];
    knots2 = [0; 0; 0; 2; 2; 2];
    
    
    %----------------------------------------------------------------------
    %  Set the initial nodes array
    %----------------------------------------------------------------------
    % Some useful constants
    plate_middle = (hole_radius + plate_length) / 2;
    constant_sq2m1 = sqrt(2) - 1;
    w23 = 1/2 + sqrt(2)/4;
    
    nodes = [0                            , hole_radius                  , 1; ...
             constant_sq2m1 * hole_radius , hole_radius                  , w23; ...
             hole_radius                  , constant_sq2m1 * hole_radius , w23; ...
             hole_radius                  , 0                            , 1; ...
             ...
             0                            , plate_middle                 , 1; ...
             constant_sq2m1 * plate_middle, plate_middle                 , w23; ...
             plate_middle                 , constant_sq2m1 * plate_middle, w23; ...
             plate_middle                 , 0                            , 1; ...
             ...
             0                            , plate_length                 , 1; ...
             plate_length                 , plate_length                 , w23; ...
             plate_length                 , plate_length                 , w23; ...
             plate_length                 , 0                            , 1];
    
    % Project the NURBS nodes to B-spline nodes
    nodes = project_up(nodes);
    
    
    %----------------------------------------------------------------------
    %  Set the initial degrees
    %----------------------------------------------------------------------
    p1 = 2;
    p2 = 2;
    
    
    %----------------------------------------------------------------------
    %  Knot refine in directions 1 and 2
    %  
    %  As the base case, we insert the knot 1 in direction 2 so that we
    %  start with 4 elements. numRefinements indicates how many more times
    %  we want to refine the knot vectors.
    %----------------------------------------------------------------------
    % Specify the knots for insertion
    temp1 = setdiff((0 : 0.5^numRefinements : 1)', [0; 1]);
    temp2 = setdiff((0 : 0.5^numRefinements : 2)', [0; 2]);
    
    knotsForInsertion1 = [(0 + temp1); (1 + temp1)];
    knotsForInsertion2 = temp2;
    
    clear temp1 temp2;
    
    
    % Refine the knot vectors
    [knots1, knots2, nodes] = ...
        refine_h_surface(knots1, knots2, nodes, p1, p2, knotsForInsertion1, knotsForInsertion2);
    
    clear knotsForInsertion1 knotsForInsertion2;
    
    
    % Project the B-spline nodes to NURBS nodes
    nodes = project_down(nodes);
    numNodes = size(nodes, 1);
    
    %{
    % Draw the NURBS patch (debugging)
    draw_nurbs_surface(knots1, knots2, nodes, p1, p2);
    axis image;
    %}
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extraction matrices
    %----------------------------------------------------------------------
    [bezierExtractions1, nodeIndexShifts1, numElements1, elementSizes1] = build_bezier_extraction(knots1, p1);
    [bezierExtractions2, nodeIndexShifts2, numElements2, elementSizes2] = build_bezier_extraction(knots2, p2);
    
    
    %----------------------------------------------------------------------
    %  Build the IEN array
    %----------------------------------------------------------------------
    IEN_array = build_ien_array(knots1, knots2, [], p1, p2, [], nodeIndexShifts1, nodeIndexShifts2, []);
    
    clear nodeIndexShifts1 nodeIndexShifts2;
    
    
    %----------------------------------------------------------------------
    %  Set the BCs_displacement array
    %  [node index, dof index, value]
    %----------------------------------------------------------------------
    % Number of nodes on an element along directions 1 and 2
    numNodesPerElement1 = p1 + 1;
    numNodesPerElement2 = p2 + 1;
    
    % Number of degrees of freedom on each element
    numDOFsPerElement = numDOFsPerNode * (numNodesPerElement1 * numNodesPerElement2);
    
    
    % Displacement of the left side of the quarter plate (fixed in x)
    value_left = 0;
    numBCs_left = (1 * numElements2) * (1 * numNodesPerElement2);
    
    % Displacement of the bottom side of the quarter plate (fixed in y)
    value_bottom = 0;
    numBCs_bottom = (1 * numElements2) * (1 * numNodesPerElement2);
    
    % Initialize the displacement BCs array
    numBCs_displacement = numBCs_left + numBCs_bottom;
    BCs_displacement = zeros(numBCs_displacement, 3);
    
    
    % Set the displacement BC for the left side
    elementIndex = find_element_index(1, (1 : numElements2)', [], numElements1, numElements2, []);
    
    nodeIndex = find_node_index(1, (1 : numNodesPerElement2)', [], numNodesPerElement1, numNodesPerElement2, [], IEN_array(:, elementIndex));
    
    BCs_displacement((1 : numBCs_left)', :) = [nodeIndex, 1 * ones(numBCs_left, 1), value_left * ones(numBCs_left, 1)];
    
    clear elementIndex nodeIndex;
    
    
    % Set the displacement BC for the bottom side
    elementIndex = find_element_index(numElements1, (1 : numElements2)', [], numElements1, numElements2, []);
    
    nodeIndex = find_node_index(numNodesPerElement1, (1 : numNodesPerElement2)', [], numNodesPerElement1, numNodesPerElement2, [], IEN_array(:, elementIndex));
    
    BCs_displacement(((numBCs_left + 1) : numBCs_displacement)', :) = [nodeIndex, 2 * ones(numBCs_bottom, 1), value_bottom * ones(numBCs_bottom, 1)];
    
    clear elementIndex nodeIndex;
    
    
    %----------------------------------------------------------------------
    %  Set the BCs_force array
    %  [node index, dof index, value]
    %----------------------------------------------------------------------
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    
    
    % Set the quadrature rule
    [z1, z2, w] = set_2d_gauss_quadrature_for_bernstein(8);
    
    % Evaluate the Bernstein polynomials and their derivatives at the
    % quadrature points for direction 1
    Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints(1));
    Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints(1));
    
    for j = 1 : numQuadraturePoints(1)
        temp = eval_1d_bernstein_der(z1(j), p1);
        
        Bernstein1_der0(:, j) = temp(:, 1);
        Bernstein1_der1(:, j) = temp(:, 2);
    end
    
    clear z1 temp;
    
    % Evaluate the Bernstein polynomials and their derivatives at the
    % quadrature points for direction 2
    Bernstein2_der0 = zeros(constant_p2p1, numQuadraturePoints(2));
    Bernstein2_der1 = zeros(constant_p2p1, numQuadraturePoints(2));
    
    for j = 1 : numQuadraturePoints(2)
        temp = eval_1d_bernstein_der(z2(j), p2);
        
        Bernstein2_der0(:, j) = temp(:, 1);
        Bernstein2_der1(:, j) = temp(:, 2);
    end
    
    clear z2 temp;
    
    
    % Tractions on the top side of the quarter plate (in x and y)
    numBCs_top = 2 * ((numElements1 / 2) * 1) * (numNodesPerElement1 * 1);
    
    % Tractions on the right side of the quarter plate (in x and y)
    numBCs_right = 2 * ((numElements1 / 2) * 1) * (numNodesPerElement1 * 1);
    
    % Initialize the force BCs array
    numBCs_force = numBCs_top + numBCs_right;
    BCs_force = zeros(numBCs_force, 3);
    
    
    % Set the force BC for the top side
    elementLocalIndex1 = (1 : (numElements1 / 2))';
    elementLocalIndex2 = numElements2;
    elementIndex = find_element_index(elementLocalIndex1, elementLocalIndex2, [], numElements1, numElements2, []);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', numNodesPerElement2, [], numNodesPerElement1, numNodesPerElement2, [], IEN_array(:, elementIndex));
    nodeIndex = kron(nodeIndex, ones(2, 1));
    
    for e2 = elementLocalIndex2
        % Find the Bezier extraction matrix
        bezierExtractions2_e = bezierExtractions2(:, :, e2);
        
        % Evaluate deta/dt (constant)
        deta_dt = elementSizes2(e2);
        dt_deta = 1 / deta_dt;
        
        % Evaluate the univariate B-splines
        Bspline2_der0 = bezierExtractions2_e * Bernstein2_der0;
        Bspline2_der1 = (dt_deta * bezierExtractions2_e) * Bernstein2_der1;
        
        
        for e1 = elementLocalIndex1
            % Find the Bezier extraction matrix
            bezierExtractions1_e = bezierExtractions1(:, :, e1);
            
            % Evaluate dxi/dt (constant)
            dxi_dt = elementSizes1(e1);
            dt_dxi = 1 / dxi_dt;
            
            % Evaluate the univariate B-splines
            Bspline1_der0 = bezierExtractions1_e * Bernstein1_der0;
            Bspline1_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
            
            
            %--------------------------------------------------------------
            %  Loop over the quadrature points
            %--------------------------------------------------------------
            % Initialize the element matrix and element RHS vector
            f_e = zeros(numDOFsPerElement, 1);
            
            % Find the positions and weights of the nodes
            nodes_e = nodes(IEN_array(:, e), [1; 2]);
            w_e = nodes(IEN_array(:, e), 3);
            
            % Evaluate the B-splines at the points
            % Matrix: (p1 + 1) * (p2 + 1) x (numQuadraturePoints1 * numQuadraturePoints2)
            Bspline_der00 = kron(Bspline2_der0, Bspline1_der0);
            Bspline_der10 = kron(Bspline2_der0, Bspline1_der1);
            Bspline_der01 = kron(Bspline2_der1, Bspline1_der0);
            
            [NURBS_der00, NURBS_der10, NURBS_der01] = eval_2d_nurbs_der1(w_e, Bspline_der00, Bspline_der10, Bspline_der01);
            
            % Evaluate the map derivatives dx/dxi at the points
            % Matrix: 2 x (numQuadraturePoints1 * numQuadraturePoints2)
            dx_dxi = nodes_e' * NURBS_der10;
            dx_deta = nodes_e' * NURBS_der01;
            
            
            for j = 1 : numQuadraturePoints(1) * numQuadraturePoints(2)
                % Form the Jacobian matrix
                JacobianMatrix = [dx_dxi(:, j)'; ...
                                  dx_deta(:, j)'];
                
                % Compute the Jacobian
                Jacobian = det(JacobianMatrix) * dxi_dt * deta_dt;
                if (Jacobian <= 0)
                   fprintf('Error: Jacobian is not positive for e = %d, j = %d\n\n', e, j); 
                end
                
                % Derivatives of NURBS with respect to the parametric
                % domain
                % Matrix: 2 x (p1 + 1) * (p2 + 1)
                NURBS_der1 = [NURBS_der10(:, j)'; ...
                              NURBS_der01(:, j)'];
            end
        end
    end
    
    
    
    % Set the force BC for the right side
    elementIndex = find_element_index(((numElements1 / 2 + 1) : numElements1)', numElements2, [], numElements1, numElements2, []);
    
    nodeIndex = find_node_index((1 : numNodesPerElement1)', numNodesPerElement2, [], numNodesPerElement1, numNodesPerElement2, [], IEN_array(:, elementIndex));
    nodeIndex = kron(nodeIndex, ones(2, 1));
    
    
    return
    
    
    %----------------------------------------------------------------------
    %  Build the ID, LM, and BC arrays
    %----------------------------------------------------------------------
    ID_array = build_id_array(BCs_displacement, numNodes, numDOFsPerNode);
    
    LM_array = build_lm_array(IEN_array, ID_array);
    
    [BCU_array, BCF_array] = build_bc_array(BCs_displacement, BCs_force, ID_array);
    
    clear BCs_displacement BCs_force;
    
    
    % Number of degrees of freedom on the patch
    numDOFs = numDOFsPerNode * numNodes;
    
    % Number of unknown degrees of freedom
    numUnknownDOFs = numDOFs - size(BCU_array, 1);
    
    
    %----------------------------------------------------------------------
    %  Store the final values for the patch
    %----------------------------------------------------------------------
    %{
    patches(1).description = description;
    patches(1).knots1 = knots1;
    patches(1).knots2 = knots2;
    patches(1).knots3 = knots3;
    patches(1).nodes = nodes;
    patches(1).p1 = p1;
    patches(1).p2 = p2;
    patches(1).p3 = p3;
    patches(1).bezierExtractions1 = bezierExtractions1;
    patches(1).bezierExtractions2 = bezierExtractions2;
    patches(1).bezierExtractions3 = bezierExtractions3;
    patches(1).numElements1 = numElements1;
    patches(1).numElements2 = numElements2;
    patches(1).numElements3 = numElements3;
    patches(1).elementSizes1 = elementSizes1;
    patches(1).elementSizes2 = elementSizes2;
    patches(1).elementSizes3 = elementSizes3;
%   patches(1).sharedNodesIndex = [];
    patches(1).numDOFsPerNode = numDOFsPerNode;
    patches(1).numDOFsPerElement = numDOFsPerElement;
    patches(1).numDOFs = numDOFs;
    patches(1).numUnknownDOFs = numUnknownDOFs;
    patches(1).IEN_array = IEN_array;
    patches(1).ID_array = ID_array;
    patches(1).LM_array = LM_array;
    patches(1).BCU_array = BCU_array;
    patches(1).BCF_array = BCF_array;
    
    clear description;
    clear knots1 knots2 knots3;
    clear nodes;
    clear p1 p2 p3;
    clear bezierExtractions1 bezierExtractions2 bezierExtractions3;
    clear numElements1 numElements2 numElements3;
    clear elementSizes1 elementSizes2 elementSizes3;
    clear numDOFsPerNode numNodesPerElement1 numNodesPerElement2 numNodesPerElement3 numDOFs numUnknownDOFs;
    clear IEN_array ID_array LM_array BCU_array BCF_array;
    %}
            
    
    %----------------------------------------------------------------------
    %  Save the assembly file
    %----------------------------------------------------------------------
    save(file_assembly, 'description', ...
                        'knots1', 'knots2', ...
                        'nodes', ...
                        'p1', 'p2', ...
                        'bezierExtractions1', 'bezierExtractions2', ...
                        'numElements1', 'numElements2', ...
                        'elementSizes1', 'elementSizes2', ...
                        'numDOFsPerNode', 'numDOFsPerElement', 'numDOFs', 'numUnknownDOFs', ...
                        'IEN_array', 'ID_array', 'LM_array', 'BCU_array', 'BCF_array', ...
                        'material_E', 'material_nu', 'material_thickness', 'numQuadraturePoints', ...
                        '-v7.3');
end


function [u_x, u_y] = find_exact_displacement(x, y, material_E, material_nu, stress_at_infinity, hole_radius)
    % Change to polar coordinates
    r = sqrt(x^2 + y^2);
    theta = atan(y / x) + pi/2;
    r_normInv = hole_radius / r;
    
    % Assume plane stress
    kappa = (3 - material_nu) / (1 + material_nu);
    mu = material_E / (2 * (1 + material_nu));
    
    % Exact displacements
    u_x =  stress_at_infinity * hole_radius / (8*mu) * (1/r_normInv * (kappa - 3) * sin(theta) + 2*r_normInv * ((-kappa + 1) * sin(theta) + sin(3*theta)) - 2*r_normInv^3 * sin(3*theta));
    u_y = -stress_at_infinity * hole_radius / (8*mu) * (1/r_normInv * (kappa + 1) * cos(theta) + 2*r_normInv * (( kappa + 1) * cos(theta) + cos(3*theta)) - 2*r_normInv^3 * cos(3*theta));
end


function [sigma_xx, sigma_yy, sigma_xy] = find_exact_stress(x, y, stress_at_infinity, hole_radius)
    % Change to polar coordinates
    r = sqrt(x^2 + y^2);
    theta = atan(y / x) + pi/2;
    r_normInv = hole_radius / r;
    
    sigma_xx =  stress_at_infinity * (  - r_normInv^2 * (0.5 * cos(2*theta) - cos(4*theta)) - 1.5*r_normInv^4 * cos(4*theta));
    sigma_yy =  stress_at_infinity * (1 - r_normInv^2 * (1.5 * cos(2*theta) + cos(4*theta)) + 1.5*r_normInv^4 * cos(4*theta));
    sigma_xy = -stress_at_infinity * (  - r_normInv^2 * (0.5 * sin(2*theta) + sin(4*theta)) + 1.5*r_normInv^4 * sin(4*theta));
end