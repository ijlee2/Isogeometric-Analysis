%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the derivatives of the 1D B-splines of degree p
%  that are nonzero at the given location.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open.
%  
%  This routine will be inefficient for finite element analysis. Please use
%  bezier_extraction and eval_bernstein_der_all to compute the derivatives
%  of the Bernstein polynomials, and find the derivatives of the B-splines
%  using the affine map between the Bernstein and B-spline domains.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      output = eval_1d_bspline_der(knots, xi, p);
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
%  1. (p + 1) x (p + 1) matrix
%  
%      The first column corresponds to the 0th derivatives of the B-splines,
%      the second column to the 1st derivatives, and so on.
%--------------------------------------------------------------------------
function output = eval_1d_bspline_der(knots, xi, p)
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
    
    % The two-dimensional array ndu contains information about the basis
    % functions evaluated at xi (stored in the upper triangle), and that
    % about the knot differences (stored in the lower triangle)
    ndu(1, 1) = 1;
    
    for j = 1 : p
        saved = 0;
        
        for k = 1 : j
            % Form the lower triangle of ndu
            ndu(j + 1, k) = right(k) + left(j - k + 1);
            temp = ndu(k, j) / ndu(j + 1, k);
            
            % Form the upper triangle of ndu
            ndu(k, j + 1) = saved + right(k) * temp;
            saved = left(j - k + 1) * temp;
        end
        
        ndu(j + 1, j + 1) = saved;
    end
    
    
    %----------------------------------------------------------------------
    %  Compute the derivatives (j is the B-spline index)
    %----------------------------------------------------------------------
    % Load the basis functions
    output = zeros(p + 1, p + 1);
    output(1 : p + 1, 1) = ndu(1 : p + 1, p + 1);
    
    for j = 0 : p
        % Alternate rows in array alpha
        s1 = 0;
        s2 = 1; 
        
        % Base case
        alpha(1, 1) = 1;
        
        % Compute the k-th derivative
        for k = 1 : p
            d = 0;
            jk = j - k;
            pk = p - k;
            
            if (jk >= 0)
                alpha(s2 + 1, 1) = alpha(s1 + 1, 1) / ndu(pk + 2, jk + 1);
                
                d = alpha(s2 + 1, 1) * ndu(jk + 1, pk + 1);
            end
            
            if (jk >= -1)
                m1 = 1;
            else
                m1 = -jk;
            end
            
            if (j - 1 <= pk)
                m2 = k - 1;
            else 
                m2 = p - j;
            end
            
            for m = m1 : m2
                alpha(s2 + 1, m + 1) = (alpha(s1 + 1, m + 1) - alpha(s1 + 1, m)) / ndu(pk + 2, jk + m + 1);
                
                d = d + alpha(s2 + 1, m + 1) * ndu(jk + m + 1, pk + 1);
            end
            
            if (j <= pk)
                alpha(s2 + 1, k + 1) = -alpha(s1 + 1, k) / ndu(pk + 2, j + 1);
                
                d = d + alpha(s2 + 1, k + 1) * ndu(j + 1, pk + 1);
            end
            
            output(j + 1, k + 1) = d;
            
            % Swap s1 and s2
            m = s1;
            s1 = s2;      
            s2 = m;
        end
    end
    
    % Multiply through by the correct factors
    coefficient = p;
    
    for k = 1 : p
        for j = 0 : p
            output(j + 1, k + 1) = coefficient * output(j + 1, k + 1);
        end
        
        coefficient = coefficient * (p - k);
    end
end
