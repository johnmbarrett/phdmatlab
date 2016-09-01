%SUMOFPROD Sum of products of all subsets
%   z = SUMOFPROD(X) where x is a column vector computes the sum of 
%   products of all subsets of the elements of x.  If X is a matrix, the
%   calculation is performed separately for each column.
%
%   [z,Y] = SUMOFPROD(X) also returns a matrix Y whose columns are the
%   coefficients of polynomials with roots specified by the columns of X.
%
%   Equivalent to sum(abs(poly(X))) when X is a vector, but seems to be
%   slightly faster when X is longer than about 16000 elements.

%   Written by John Barrett (john.barrett@cantab.net)
%   $Revision: 1.0 $  $Date: 2013/11/29 17:26:50 $
function [z,Y] = sumofprod(X)
    Y = zeros(size(X)+[1 0]);
    z = zeros(1,size(X,2));
    
    for ii = 1:size(X,2)
        x = X(:,ii);
        y = 1;
        
        % By Vieta's formulas (http://en.wikipedia.org/wiki/Vieta%27s_formulas),
        % if the polynomial (x-x1)(x-x2)...(x-xn) with roots x1, x2, ..., xn
        % has coefficients an, ..., a0, then the (n-k)th coefficient is the sum
        % of the products of all k-element subsets of the roots (an = 1 = sum
        % of the empty product)
        for jj = 1:numel(x)
            % convolution is equivalent to polynomial multiplication (see
            % Matlab documentation for conv)
            y = conv(y,[1 -x(jj)]);
        end
    
        Y(:,ii) = y;
        z(ii) = sum(abs(y));
    end
end