classdef BSpline
    %BSPLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degree
        pieces
    end
    
    methods
        function obj = BSpline(degree)
            obj.degree = degree;
            cur(degree+1) = Polynomial(0.0);
            cur(1) = Polynomial(1.0);
            for k=1:degree
                prev = cur;
                kr = 1 / k;
                cur(1) = Polynomial([0;kr]) * prev(1);
                for i=1:k-1
                    cur(i+1) = Polynomial([i*kr;kr]) * prev(i+1) + Polynomial([1-(i-1)*kr; -kr]) * prev(i);
                end
                cur(k+1) = Polynomial([kr;-kr]) * prev(k);
            end
            obj.pieces = cur;
        end
    end
end

