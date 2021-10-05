classdef PPHelper < handle
    methods(Static)
        % Defines a piecewise linear function such that x_i |-> y_i for
        % y_i scalar
        function pp = makePiecewiseLinear(x,y)
            y = y(:);
            x = x(:);
            y_const = y(1:(end-1));
            y_lin = (y(2:end) - y(1:(end-1))) ./ (x(2:end )- x(1:(end-1)));
            pp = mkpp(x, [y_lin y_const]);
        end
        
        % Defines a piecewise linear function such that x_i |-> y_i for
        % y_i = y(:,i) a vector.
        function pp = makePiecewiseLinearVector(x,y)
            x = x(:)';
            y_const = y(:,1:(end-1));
            y_lin = (y(:,2:end) - y(:,1:(end-1))) ./ (x(2:end) - x(1:(end-1)));
            pp = mkpp(x, cat(3, y_lin, y_const), size(y,1));
        end
        
        % Defines a piecewise linear function such that x_i |-> y_i for
        % y_i = y(:,:,i) a matrix.
        function pp = makePiecewiseLinearMatrix(x,y)
            x = x(:);
            y_const = y(:,:,1:(end-1));
            y_lin = (y(:,:,2:end) - y(:,:,1:(end-1))) ./ permute(x(2:end )- x(1:(end-1)), [2 3 1]);
            pp = mkpp(x, cat(4, y_lin, y_const), [size(y,1) size(y,2)]);
        end
    end
end