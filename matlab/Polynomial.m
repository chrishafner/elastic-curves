classdef Polynomial < handle
    properties
        coeff
    end
    
    methods
        function obj = Polynomial(coeff)
            if nargin < 1
                coeff = [];
            end
            obj.coeff = coeff(:);
            obj.truncateZeros();
        end
        
        function truncateZeros(obj)
            n = find(obj.coeff ~= 0, 1, 'last');
            if isempty(n)
                obj.coeff = 0;
            else
                obj.coeff = obj.coeff(1:n);
            end
        end
        
        function ret = evaluate(obj, t)
            szt = size(t);
            t = t(:)';
            exps = (0:length(obj.coeff)-1)';
            ret = reshape(sum(obj.coeff .* (t.^exps), 1), szt);
        end
        
        function ret = evaluateD(obj, t)
            szt = size(t);
            t = t(:)';
            exps = (0:length(obj.coeff)-2)';
            ret = reshape(sum((obj.coeff(2:end) .* (exps+1)) .* (t.^exps), 1), szt);
        end
        
        function ret = evaluateD2(obj, t)
            szt = size(t);
            t = t(:)';
            exps = (0:length(obj.coeff)-3)';
            fac = (exps+1).*(exps+2);
            ret = reshape(sum((obj.coeff(3:end) .* fac) .* (t.^exps), 1), szt);
        end
        
        function ret = evaluateD3(obj, t)
            szt = size(t);
            t = t(:)';
            exps = (0:length(obj.coeff)-4)';
            fac = (exps+1).*(exps+2).*(exps+3);
            ret = reshape(sum((obj.coeff(4:end) .* fac) .* (t.^exps), 1), szt);
        end
        
        function ret = mtimes(a,b)
            if isnumeric(a) && isa(b, 'Polynomial')
                ret = Polynomial(b.coeff * a);
            elseif isnumeric(b) && isa(a, 'Polynomial')
                ret = Polynomial(a.coeff * b);
            elseif isa(a, 'Polynomial') && isa(b, 'Polynomial')
                an = length(a.coeff);
                bn = length(b.coeff);
                d = an + bn - 1;
                c = zeros(d, 1);
                for i=1:an
                    c(i:i+bn-1) = c(i:i+bn-1) + a.coeff(i) * b.coeff;
                end
                ret = Polynomial(c);
            end
        end
        
        function ret = mrdivide(a,b)
            if isnumeric(b) && isa(a, 'Polynomial')
                ret = Polynomial(a.coeff / b);
            end
        end
        
        function ret = plus(a,b)
            if isnumeric(a) && isa(b, 'Polynomial')
                ret = Polynomial([b.coeff(1) + a; b.coeff(2:end)]);
            elseif isnumeric(b) && isa(a, 'Polynomial')
                ret = Polynomial([a.coeff(1) + b; a.coeff(2:end)]);
            elseif isa(a, 'Polynomial') && isa(b, 'Polynomial')
                an = length(a.coeff);
                bn = length(b.coeff);
                if an > bn
                    c = a.coeff;
                    c(1:bn) = c(1:bn) + b.coeff;
                else
                    c = b.coeff;
                    c(1:an) = c(1:an) + a.coeff;
                end
                ret = Polynomial(c);
            end
        end
        
        function ret = minus(a,b)
            ret = a + (-b);
        end
        
        function ret = uminus(a)
            ret = Polynomial(-a.coeff);
        end
        
        function ret = uplus(a)
            ret = Polynomial(a.coeff);
        end
    end
end