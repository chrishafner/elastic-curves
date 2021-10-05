classdef NumOpt < handle
    methods(Static)
        function [ret, converged] = newtonsMethod(f_fun, df_fun, x, threshold, max_iter)
            f = f_fun(x);
            iter = 0;
            while (abs(f) > threshold && iter < max_iter)
                df = df_fun(x);
                x = x - f/df;
                f = f_fun(x);
                iter = iter + 1;
            end
            ret = x;
            converged = abs(f) <= threshold;
        end
    end
end