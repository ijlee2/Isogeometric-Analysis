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
%      initialize(path_to_assembly_directory, path_to_results_directory)
%  
%  where,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement and phase fields (.mat files)
%--------------------------------------------------------------------------
function initialize(path_to_assembly_directory, path_to_results_directory)
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'              , ...
         'numDOFs'                 , ...
         ...
         'p_exp'                   , ...
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
            u1(LM_array1(offset + i, e)) = function_u_exact(x(i), u_L);
            u2(LM_array2(offset + i, e)) = function_c_exact(x(i), material_ell_0, p_exp);
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
    save(sprintf('%sfile_results_alternation%d', path_to_results_directory, 0), 'u1', 'u2', '-v6');
end


function u = function_u_exact(x, u_L)
    if (x < 0)
        u = 0;
    else
        u = u_L;
    end
end


function c = function_c_exact(x, material_ell_0, p_exp)
    % Normalize the physical domain
    y = abs(x)/material_ell_0;
    
    c = 1 - exp(-y^p_exp);
end
