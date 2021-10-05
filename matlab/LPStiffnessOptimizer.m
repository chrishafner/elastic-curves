classdef LPStiffnessOptimizer < handle
    properties
        gamma = zeros(2,0) % curve samples
        kappa = zeros(1,0) % corresponding curvature samples
        e_gravity = [0;0] % gravity vector e
        gamma_infl = zeros(2,0) % list of curve inflection points
        
        v_weights % vertex weights for integrating vertex-based quantities like stiffness
    end
    
    methods
        function obj = LPStiffnessOptimizer(varargin)
            if nargin<2
                error('LPStiffnessOptimizer requires at least two input arguments.');
            end
            
            if isnumeric(varargin{1})
                obj.gamma = varargin{1};
                obj.kappa = varargin{2};
            else
                spline = varargin{1};
                ns = varargin{2};
                
                [~, p_infl] = spline.findInflectionPoints();
                t = linspace(0, spline.t_max, ns);
                gamma = spline.evaluate(t);
                curv = spline.curvature(t);
                gamma_to = gamma(:,2:end) - gamma(:,1:end-1);
                seg_lengths = sqrt(sum(gamma_to.^2, 1));
                eff_lengths = 0.5 * [seg_lengths(1), seg_lengths(1:end-1) + seg_lengths(2:end), seg_lengths(end)];

                obj.gamma = gamma;
                obj.kappa = curv;
                obj.gamma_infl = p_infl;
                obj.v_weights = eff_lengths;
            end
        end
        
        % Optimize for stiffnesses with no inflections, and no gravity
        function [K,a,b] = optimizeSimple(obj)
            n = length(obj.kappa);
            
            K_mat = [obj.gamma ./ obj.kappa; 1./ obj.kappa]';
            A = [K_mat, zeros(n,1); K_mat, -ones(n,1)];
            sense = [repmat('>', n, 1); repmat('<', n, 1)];
            rhs = [ones(n,1); zeros(n,1)];
            objective = [0;0;0;1];
            lb = [-inf(3,1); 1];
            
            % K*kappa = a + <b,gamma>
            % Vars: b1,b2,a,M
            model = struct('A',sparse(A),'sense',sense,'rhs',rhs,'lb',lb,'obj',objective);
            params = struct('OutputFlag', 0);
            sln = gurobi(model, params);
            
            b = sln.x([1;2]);
            a = sln.x(3);
            % M = sln.x(4);
            
            K = (sum(obj.gamma .* b, 1) + a) ./ obj.kappa;
            K = K';
        end
        
        function [K,a,b,max_err] = optimizeWithTarget(obj, K0)
            n = length(obj.kappa);
            K_mat = [obj.gamma ./ obj.kappa; 1./ obj.kappa]';
            A = [K_mat, zeros(n,1); K_mat, -ones(n,1)];
            sense = [repmat('>', n, 1); repmat('<', n, 1)];
            rhs = [K0(:); K0(:)];
            objective = [0;0;0;1];
            lb = -inf(4,1);
            
            model = struct('A',sparse(A),'sense',sense,'rhs',rhs,'lb',lb,'obj',objective);
            params = struct('OutputFlag', 0);
            sln = gurobi(model, params);
            
            b = sln.x([1;2]);
            a = sln.x(3);
            max_err = sln.x(4);
            
            K = (sum(obj.gamma .* b, 1) + a) ./ obj.kappa;
            K = K';
        end
        
        function [K, a, b] = optimizeWithInflections(obj)
            n = length(obj.kappa);
            n_infl = size(obj.gamma_infl, 2);
            
            K_mat = [obj.gamma ./ obj.kappa; 1./ obj.kappa]';
            A = [K_mat, zeros(n,1); K_mat, -ones(n,1); obj.gamma_infl', ones(n_infl, 1), zeros(n_infl, 1)];
            sense = [repmat('>', n, 1); repmat('<', n, 1); repmat('=', n_infl, 1)];
            rhs = [ones(n,1); zeros(n,1); zeros(n_infl,1)];
            objective = [0;0;0;1];
            lb = [-inf(3,1); 1];
            
            % K*kappa = a + <b,gamma>
            % Vars: b1,b2,a,M
            model = struct('A',sparse(A),'sense',sense,'rhs',rhs,'lb',lb,'obj',objective);
            params = struct('OutputFlag', 0);
            sln = gurobi(model, params);
            
            if strcmp(sln.status, 'OPTIMAL')
                b = sln.x([1;2]);
                a = sln.x(3);
                % M = sln.x(4);

                K = (sum(obj.gamma .* b, 1) + a) ./ obj.kappa;
                K = K';
            else
                K = [];
                a = [];
                b = [];
            end
        end
        
        function [K, a, b] = optimizeWithGravity(obj)
            n = length(obj.kappa);
            
            % Indices of variables (column indices in A)
            ia = 1;
            ib = [2;3];
            iK = (4:n+3)';
            iKint = (n+4:2*n+3)';
            iP = (2*n+4:3*n+3)'; % P = int (K <gamma, R^T e>)
            iM = 3*n+4;
            
            nv = 3*n + 4;
            lb = -inf(nv,1);
            lb(iK) = 1;
            
            objective = zeros(nv,1);
            objective(iM) = 1;
            
            % Equations and inequalities (number of eq.)
            % 1) K <= M (n)
            % 2) Definition of Kint (n)
            % 3) Definition of P (n)
            % 4) Moment equilibrium (n)
            
            roff = 0; % row offset (row offset for (in)equalities in A)
            
            % K<=M -> K-M <= 0
            ci1 = [iK, iM*ones(n,1)];
            ri1 = repmat((roff+1:roff+n)', 1, 2);
            val1 = [ones(n,1), -ones(n,1)];
            sense1 = repmat('<',n,1);
            rhs1 = zeros(n,1);
            
            roff = roff + n;
            
            % Kint(i) = Kint(i+1) + K(i)*v_weights(i)
            % Kint(n) = 0 + K(n)*v_weights(n)
            ci2 = [iKint, [iKint(2:end); iKint(end)], iK];
            ri2 = repmat((roff+1:roff+n)', 1, 3);
            val2 = [ones(n,1), [-ones(n-1,1); 0], -obj.v_weights'];
            sense2 = repmat('=',n,1);
            rhs2 = zeros(n,1);
            
            roff = roff + n;
            
            % P(i) = P(i+1) + K(i)*<gamma,R^T e>*v_weights(i)
            % P(end) = 0 + K(n)*<gamma,R^T e>*v_weights(n)
            Rte = [obj.e_gravity(2); -obj.e_gravity(1)];
            f = sum(obj.gamma .* Rte, 1);
            
            ci3 = [iP, [iP(2:end); iP(end)], iK];
            ri3 = repmat((roff+1:roff+n)', 1, 3);
            val3 = [ones(n,1), [-ones(n-1,1); 0], -f' .* obj.v_weights'];
            sense3 = repmat('=',n,1);
            rhs3 = zeros(n,1);
            
            roff = roff + n;
            
            % -K*kappa + <gamma,b> + f*Kint - P + a = 0
            ci4 = [iK, repmat(ib',n,1), iKint, iP, repmat(ia,n,1)];
            ri4 = repmat((roff+1:roff+n)', 1, 6);
            val4 = [-obj.kappa', obj.gamma', f', -ones(n,1), ones(n,1)];
            sense4 = repmat('=',n,1);
            rhs4 = zeros(n,1);
            
            roff = roff + n;
            
            
            ci = [ci1(:); ci2(:); ci3(:); ci4(:)];
            ri = [ri1(:); ri2(:); ri3(:); ri4(:)];
            val = [val1(:); val2(:); val3(:); val4(:)];
            A = sparse(ri, ci, val, roff, nv);
            
            sense = [sense1(:); sense2(:); sense3(:); sense4(:)];
            rhs = [rhs1(:); rhs2(:); rhs3(:); rhs4(:)];
            
            model = struct('A',A,'sense',sense,'rhs',rhs,'lb',lb,'obj',objective);
            params = struct('OutputFlag', 0);
            sln = gurobi(model, params);
            
            K = sln.x(iK);
            a = sln.x(ia);
            b = sln.x(ib);
            
        end
        
        function [K, a, b] = fineTuneWithGravity(obj, M_opt, M_eps)
            n = length(obj.kappa);
            
            % Indices of variables (column indices in A)
            ia = 1;
            ib = [2;3];
            iK = (4:n+3)';
            iKint = (n+4:2*n+3)';
            iP = (2*n+4:3*n+3)'; % P = int (K <gamma, R^T e>)
            iM = 3*n+4;
            iKabs = (3*n+5:4*n+2)';
            
            nv = 4*n + 2;
            lb = -inf(nv,1);
            lb(iK) = 1;
            
            objective = zeros(nv,1);
%             objective(iM) = 1;
            objective(iKabs) = 1;
            
            % Equations and inequalities (number of eq.)
            % 1) K <= M (n)
            % 2) Definition of Kint (n)
            % 3) Definition of P (n)
            % 4) Moment equilibrium (n)
            % 5) M_opt*(1-M_eps) <= M <= M_opt(1+M_eps)
            % 6) iKabs >= dif(iK), iKabs >= -dif(iK)
            
            roff = 0; % row offset (row offset for (in)equalities in A)
            
            % K<=M -> K-M <= 0
            ci1 = [iK, iM*ones(n,1)];
            ri1 = repmat((roff+1:roff+n)', 1, 2);
            val1 = [ones(n,1), -ones(n,1)];
            sense1 = repmat('<',n,1);
            rhs1 = zeros(n,1);
            
            roff = roff + n;
            
            % Kint(i) = Kint(i+1) + K(i)*v_weights(i)
            % Kint(n) = 0 + K(n)*v_weights(n)
            ci2 = [iKint, [iKint(2:end); iKint(end)], iK];
            ri2 = repmat((roff+1:roff+n)', 1, 3);
            val2 = [ones(n,1), [-ones(n-1,1); 0], -obj.v_weights'];
            sense2 = repmat('=',n,1);
            rhs2 = zeros(n,1);
            
            roff = roff + n;
            
            % P(i) = P(i+1) + K(i)*<gamma,R^T e>*v_weights(i)
            % P(end) = 0 + K(n)*<gamma,R^T e>*v_weights(n)
            Rte = [obj.e_gravity(2); -obj.e_gravity(1)];
            f = sum(obj.gamma .* Rte, 1);
            
            ci3 = [iP, [iP(2:end); iP(end)], iK];
            ri3 = repmat((roff+1:roff+n)', 1, 3);
            val3 = [ones(n,1), [-ones(n-1,1); 0], -f' .* obj.v_weights'];
            sense3 = repmat('=',n,1);
            rhs3 = zeros(n,1);
            
            roff = roff + n;
            
            % -K*kappa + <gamma,b> + f*Kint - P + a = 0
            ci4 = [iK, repmat(ib',n,1), iKint, iP, repmat(ia,n,1)];
            ri4 = repmat((roff+1:roff+n)', 1, 6);
            val4 = [-obj.kappa', obj.gamma', f', -ones(n,1), ones(n,1)];
            sense4 = repmat('=',n,1);
            rhs4 = zeros(n,1);
            
            roff = roff + n;
            
            %  M >= M_opt*(1-M_eps)
            %  M <= M_opt*(1+M_eps)
            ci5 = [iM;iM];
            ri5 = [roff+1;roff+2];
            val5 = [1;1];
            sense5 = ['>'; '<'];
            rhs5 = [M_opt*(1-M_eps); M_opt*(1+M_eps)];
            
            roff = roff + 2;
            
            
            % Kabs(i) + K(i-1) - 2*K(i) + K(i+1) >= 0
            % Kabs(i) - K(i-1) + 2*K(i) - K(i+1) >= 0
            ci6 = repmat([iKabs iK(1:end-2) iK(2:end-1) iK(3:end)],2,1);
            ri6 = repmat((roff+1:roff+2*(n-2))', 1, 4);
            val6 = [repmat([1 1 -2 1],n-2,1); repmat([1 -1 2 -1],n-2,1)];
            sense6 = repmat('>', 2*(n-2), 1);
            rhs6 = zeros(2*(n-2), 1);
            
            roff = roff + 2*(n-2);
            
            ci = [ci1(:); ci2(:); ci3(:); ci4(:); ci5(:); ci6(:)];
            ri = [ri1(:); ri2(:); ri3(:); ri4(:); ri5(:); ri6(:)];
            val = [val1(:); val2(:); val3(:); val4(:); val5(:); val6(:)];
            A = sparse(ri, ci, val, roff, nv);
            
            sense = [sense1(:); sense2(:); sense3(:); sense4(:); sense5(:); sense6(:)];
            rhs = [rhs1(:); rhs2(:); rhs3(:); rhs4(:); rhs5(:); rhs6(:)];
            
            model = struct('A',A,'sense',sense,'rhs',rhs,'lb',lb,'obj',objective);
            params = struct('OutputFlag', 0);
            sln = gurobi(model, params);
            
            K = sln.x(iK);
            a = sln.x(ia);
            b = sln.x(ib);
            
        end
    end
end