classdef PlaneCurveSolver < handle
    properties
        position_samples
        curvature_samples
        tangent_samples
        parameter_samples
        inflection_points
        
        soln
    end
    
    methods
        function obj = PlaneCurveSolver()
            n = 100;
            ang = linspace(0,pi,n);
            b = 0.5;
            obj.position_samples = [b * cos(ang);sin(ang)];
            to = obj.position_samples(:,2:end) - obj.position_samples(:,1:end-1);
            lens = sqrt(sum(to.^2,1));
            obj.parameter_samples = [0 cumsum(lens)];
            obj.tangent_samples = [-b * sin(ang);cos(ang)];
            obj.tangent_samples = obj.tangent_samples ./ sqrt(sum(obj.tangent_samples.^2,1));
            obj.curvature_samples = b ./ ((b*sin(ang).^2 + cos(ang).^2).^(3/2));
            %obj.curvature_samples = ones(20,1);
            obj.inflection_points = zeros(2,0);
        end
        
        function computeStiffness(obj,a,b)
            n = size(obj.position_samples, 2);
            K = zeros(n,1);
            for i=1:n
                curv = obj.curvature_samples(i);
                K(i) = (a + dot(b, obj.position_samples(:,i))) / curv;
            end
            if min(K) <= 0
                if max(K) >= 0
                    fprintf('Stiffness infeasible!');
                else
                    K = -K;
                    a = -a;
                    b = -b;
                end
            end
            obj.soln = struct('K',K,'a',a,'b',b, ...
                'position_samples',obj.position_samples,'curvature_samples',obj.curvature_samples,...
                'inflection_points',obj.inflection_points,'tangent_samples',obj.tangent_samples);
        end
        
        function solve(obj)
            n = size(obj.position_samples, 2);
            n_infl = size(obj.inflection_points, 2);
            Ki = (1:n)';
            ai = n+1;
            bi = [n+2; n+3];
            Mi = n+4;
            
            rowi = zeros(4*n + 3*n_infl + 2*n, 1);
            coli = zeros(size(rowi));
            values = zeros(size(rowi));
            for i=1:n
                coli((4*i-3):(4*i)) = [Ki(i); ai; bi];
                rowi((4*i-3):(4*i)) = i;
                values((4*i-3):(4*i)) = [obj.curvature_samples(i); -1; -obj.position_samples(:,i)];
            end
            offset = 4*n;
            row_offset = n;
            for i=1:n_infl
                coli(offset + ((3*i-2):(3*i))) = [ai; bi];
                rowi(offset + ((3*i-2):(3*i))) = row_offset + i;
                values(offset + ((3*i-2):(3*i))) = [1; obj.inflection_points(:,i)]; 
            end
            offset = offset + 3*n_infl;
            row_offset = row_offset + n_infl;
            for i=1:n
                coli(offset + ((2*i-1):(2*i))) = [Ki(i); Mi];
                rowi(offset + ((2*i-1):(2*i))) = row_offset + i;
                values(offset + ((2*i-1):(2*i))) = [1; -1];
            end
            
            sense = char(2*n+n_infl, 1);
            sense(1:n) = '=';
            sense((n+1):(n+n_infl)) = '=';
            sense((n+n_infl+1):(n+n_infl+n)) = '<';
            
            lb = -inf(n+4, 1);
            lb(Ki) = 1.;
            
            fobj = zeros(n+4, 1);
            fobj(Mi) = 1.;
            
            A = sparse(rowi, coli, values, 2*n+n_infl, n+4);
            model = struct('A', A, 'obj', fobj, 'lb', lb, 'sense', sense);
            
            result = gurobi(model);
            K = result.x(Ki);
            a = result.x(ai);
            b = result.x(bi);
            
            obj.soln = struct('K',K,'a',a,'b',b, ...
                'position_samples',obj.position_samples,'curvature_samples',obj.curvature_samples,...
                'inflection_points',obj.inflection_points,'tangent_samples',obj.tangent_samples);
        end
    end
end