%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine initializes the coefficients for the displacement and
%  phase field. The coefficients are set to be the values of the exact
%  solutions at the nodes.
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
         'numDOFs'                 , ...
         ...
         'material_ell_0');
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'elementType'       , ...
         'p1'                , ...
         'nodes'             , ...
         'numElements1'      , ...
         'IEN_array');
    
    % Load the BCs
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'LM_array1' , 'LM_array2' , ...
         'BCU_array1', 'BCU_array2', ...
         'u_L');
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over elements (e = e1)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constant
    constant_p1p1 = p1 + 1;
    
    % Initialize the solution vectors
    u1 = zeros(numDOFs, 1);
    u2 = zeros(numDOFs, 1);
    
    % Index that was last used to set the entry of u
    lastIndex_for_u = 0;
    
    
    for e = 1 : numElements1
        %------------------------------------------------------------------
        %  Evaluate the exact solutions
        %------------------------------------------------------------------
        % Find the positions of the nodes
        if (e == 1)
            x = nodes(IEN_array(:, e), :);
        else
            x = setdiff(nodes(IEN_array(:, e), :), nodes_em1);
        end
        
        
        numNodesToConsider = size(x, 1);
        offset = constant_p1p1 - numNodesToConsider;
        
        for i = 1 : numNodesToConsider
            u1(LM_array1(offset + i, e)) = function_u_exact(x(i), u_L, material_ell_0, order);
            u2(LM_array2(offset + i, e)) = function_c_exact(x(i), material_ell_0, order);
        end
        
        
        % Save the nodes on this element for comparison against the nodes
        % on the next element
        nodes_em1 = nodes(IEN_array(:, e));
        
        % Update the index
        lastIndex_for_u = lastIndex_for_u + numNodesToConsider;
    end
    
    
    % Check for error
    if (lastIndex_for_u ~= numDOFs)
        fprintf('\n');
        fprintf('ERROR: The DOFs were not initialized correctly.\n\n');
        
        quit;
    end
    
    % Prescribe the known coefficients
    u1(BCU_array1(:, 1)) = BCU_array1(:, 2);
    u2(BCU_array2(:, 1)) = BCU_array2(:, 2);
    
    
    %----------------------------------------------------------------------
    %  Save the fields
    %----------------------------------------------------------------------
    switch order
        case {2, 4}
            u2 = [u2; 0];
            
        case {6, 8}
            u2 = [u2; 0; 0];
            
    end
    
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
