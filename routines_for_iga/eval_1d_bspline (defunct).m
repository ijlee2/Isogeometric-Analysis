%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  This code is based on that by Piegl and Tiller in their NURBS book.
%  
%  
%  Summary:
%  
%  This routine evaluates the 1D B-splines of degree p that are nonzero
%  at the given location.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open.
%  
%  This routine will be inefficient for finite element analysis. Please use
%  bezier_extraction and eval_bernstein_all to compute the Bernstein
%  polynomials, and find the B-splines using the affine map between the
%  the Bernstein and B-spline domains.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      output = eval_1d_bspline(knots, xi, p);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      xi is a knot value
%      p is the degree of the B-splines
%  
%  
%  Output:
%  
%  1. (p + 1) x 1 column vector
%  
%      The i-th row returns the evaluation of the local i-th B-spline
%      at the location t.
%--------------------------------------------------------------------------
function output = eval_1d_bspline(knots, xi, p)
    %----------------------------------------------------------------------
    %  Find the last B-spline that is nonzero at xi
    %----------------------------------------------------------------------
    numKnots = size(knots, 1);
    
    % Check that xi is indeed in the range of the knot vector Xi
    if (xi < knots(1) || xi > knots(numKnots))
        fprintf('Error: the given knot xi is outside of the range of the knot vector Xi.\n');
        return;
        
    % Special case: xi is the last the knot value in Xi
    elseif (xi == knots(numKnots))
        for lastIndex = numKnots : -1 : 1
            if (xi > knots(lastIndex))
                break;
            end
        end
        
    % Do a binary search for xi
    else
        low = 1;
        high = numKnots;
        mid = floor((low + high)/2);
        
        while (xi < knots(mid) || xi >= knots(mid + 1))
            if (xi < knots(mid))
                high = mid;
            else
                low = mid;
            end
            
            mid = floor((low + high)/2);
        end
        
        lastIndex = mid;
    end
    
    
    %----------------------------------------------------------------------
    %  Form the knot difference arrays left and right, which we define as,
    %  
    %      left(j)  = xi - Xi(lastIndex + 1 - j)
    %      right(j) = Xi(lastIndex + j) - xi
    %----------------------------------------------------------------------
    left  = xi - knots(lastIndex + 1 - (1 : p));
    right = knots(lastIndex + (1 : p)) - xi;
    
    % Base case
    output = zeros(p + 1, 1);
    output(1) = 1;
    for j = 1 : p
        saved = 0;
        
        for k = 1 : j
            temp = output(k) / (right(k) + left(j - k + 1));
            
            output(k) = saved + right(k) * temp;
            
            saved = left(j - k + 1) * temp;
        end
        
        output(j + 1) = saved;
    end
end