classdef AbsoluteAngleElastica < matlab.mixin.Copyable
    properties
        seg_lengths
        stiffnesses
        kappa_natural
        alpha_start
        alpha_end
        num_segments
        p_start
        
        total_length
        eff_lengths
        K_times_l
        K_over_l
        e_gravity = [0;0] % 2-by-1 vector, for computing gravitational energy as integral(K*<gamma,e_gravity>)
        
        % Positional constraint: Let p be the point with index pc_index.
        % Then, <pc_coeff, p> = pc_value. (This is per column of pc_index,
        % pc_coeff, and pc_value.)
        pc_index
        pc_coeff
        pc_value
        
        right_end_free = false
    end
    
    methods(Static)
        function [obj, alpha_init] = import(filename, keep_num_samples)
            file = fopen(filename);
            n = textscan(file, '%u', 1);
            n = double(n{1});
            
            if nargin < 2
                keep_num_samples = n;
            end

            seg_lengths = textscan(file, '%f', n);
            seg_lengths = seg_lengths{1}(1:keep_num_samples);
            stiffnesses = textscan(file, '%f', n+1);
            stiffnesses = stiffnesses{1}(1:keep_num_samples+1);
            p_start = textscan(file, '%f', 2);
            p_start = p_start{1};
            p_end = textscan(file, '%f', 2);
            p_end = p_end{1};
            alpha_init = textscan(file, '%f', n+2);
            alpha_init = alpha_init{1}(1:keep_num_samples+2);
            
            if keep_num_samples < n
                p_start = [0;0];
                p_end = sum([cos(alpha_init(2:end-1)), sin(alpha_init(2:end-1))] .* seg_lengths, 1)';
            end
            
            angle = atan2(p_end(2)-p_start(2), p_end(1)-p_start(1));
            magn = norm(p_end-p_start);
            p_start = [0;0];
            p_end = [magn;0];
            alpha_init = alpha_init - angle;
            
            alpha_start = alpha_init(1);
            alpha_end = alpha_init(end);
            alpha_init = alpha_init(2:end-1);
            
            obj = AbsoluteAngleElastica(seg_lengths, stiffnesses, alpha_start, alpha_end, p_start);
            obj.addEndPointConstraint(p_end);
            
            fclose(file);
        end
        
        
        function K_fixed = fixStiffnessForEquilibrium(K0, alpha, seg_lengths, alpha_start, alpha_end)
            eff_lengths = 0.5 * ...
                [seg_lengths(1);
                seg_lengths(1:end-1) + seg_lengths(2:end);
                seg_lengths(end)];
            alpha_dif = [alpha; alpha_end] - [alpha_start; alpha];
            kappa_eff = alpha_dif ./ eff_lengths;
            lsa = seg_lengths .* sin(alpha);
            lca = seg_lengths .* cos(alpha);
            
            n = length(alpha);
            nvars = n+3;
            A_vals = [kappa_eff(1:end-1); -kappa_eff(2:end); -lsa; lca];
            A_rows = repmat((1:n)', 4, 1);
            A_cols = [3:n+2, 4:n+3, ones(1,n), 2*ones(1,n)]';
            Aeq = sparse(A_rows, A_cols, A_vals, n, nvars);
            beq = zeros(n,1);
            
            C_vals = ones(n+1,1);
            C_rows = (1:n+1)';
            C_cols = (3:n+3)';
            C = sparse(C_rows, C_cols, C_vals, n+1, nvars);
            
            options = optimoptions('lsqlin', ...
                'Display', 'off');
            ls_result = lsqlin(C,K0,[],[],Aeq,beq,[],[],[],options);
            K_fixed = ls_result(3:end);
            % K_fixed = K_fixed / min(K_fixed);
        end
    end
    
    methods
        function obj = AbsoluteAngleElastica(seg_lengths, stiffnesses, alpha_start, alpha_end, p_start)
            % Elastica with n segments
            % seg_lengths: length of all n segments
            % stiffnesses: stiffnesses at all n+1 vertices
            % p_start, p_end: end points of elastica
            obj.seg_lengths = seg_lengths;
            obj.stiffnesses = stiffnesses;
            obj.kappa_natural = zeros(size(stiffnesses));
            obj.alpha_start = alpha_start;
            obj.alpha_end = alpha_end;
            obj.p_start = p_start;
            
            obj.recomputeDerived();
            
            obj.pc_index = zeros(1,0);
            obj.pc_coeff = zeros(2,0);
            obj.pc_value = zeros(1,0);
            
%             obj.p_end = p_end;
%             obj.p_to = p_end - p_start;
        end
        
        function resetConstraints(obj)
            obj.pc_index = zeros(1,0);
            obj.pc_coeff = zeros(2,0);
            obj.pc_value = zeros(1,0);
        end
        
        function recomputeDerived(obj)
            obj.num_segments = length(obj.seg_lengths);
            obj.eff_lengths = 0.5 * ...
                [obj.seg_lengths(1);
                obj.seg_lengths(1:end-1) + obj.seg_lengths(2:end);
                obj.seg_lengths(end)];
            obj.K_times_l = obj.stiffnesses .* obj.eff_lengths;
            obj.K_over_l = obj.stiffnesses ./ obj.eff_lengths;
            obj.total_length = sum(obj.seg_lengths);
        end
        
        function scaleToLength(obj, l_new)
            l = sum(obj.seg_lengths);
            factor = l_new / l;
            obj.seg_lengths = obj.seg_lengths * factor;
            obj.p_start = obj.p_start * factor;
            obj.pc_value = obj.pc_value * factor;
            obj.recomputeDerived();
        end
        
        function alpha_resampled = resampleLike(obj, other, alpha)
            s = cumsum([0; obj.seg_lengths]);
            kpp = mkpp(s,[(obj.stiffnesses(2:end)-obj.stiffnesses(1:end-1))./obj.seg_lengths, obj.stiffnesses(1:end-1)]);
            s_other = cumsum([0; other.seg_lengths]);
            k_new = ppval(kpp, s_other);
            
            if nargin > 2
                s_mid = 0.5 * (s(1:end-1) + s(2:end));
                app = mkpp(s_mid,[(alpha(2:end)-alpha(1:end-1))./obj.eff_lengths(2:end-1), alpha(1:end-1)]);
                s_other_mid = 0.5 * (s_other(1:end-1) + s_other(2:end));
                alpha_resampled = ppval(app, s_other_mid);
            end
            
            obj.pc_index = other.pc_index;
            obj.seg_lengths = other.seg_lengths;
            obj.stiffnesses = k_new;
            obj.recomputeDerived();
        end
        
        function s_mid = sMid(obj)
            s = cumsum([0; obj.seg_lengths]);
            s_mid = 0.5 * (s(1:end-1) + s(2:end));
        end
        
        function setStiffness(obj, stiffness)
            obj.stiffnesses = stiffness;
            obj.recomputeDerived();
        end
        
        function addLinearPointConstraint(obj, index, coeff, value)
            obj.pc_index = [obj.pc_index index];
            obj.pc_coeff = [obj.pc_coeff coeff(:)];
            obj.pc_value = [obj.pc_value value];
        end
        
        function addEndPointConstraint(obj, p_end)
            n = obj.num_segments;
            obj.pc_index = [obj.pc_index n n];
            obj.pc_coeff = [obj.pc_coeff eye(2)];
            obj.pc_value = [obj.pc_value p_end(:)'];
        end
        
        function ret = computeCurvature(obj, alpha)
            ret = ([alpha; obj.alpha_end] - [obj.alpha_start; alpha]) ./ obj.eff_lengths;
        end
        
        % Computes the energy, the gradient of the energy wrt all unknowns
        % (excluding Lagrangian multipliers), and the Hessian of the energy
        % in sparse-constructor form
        function [E,grad,H_row,H_col,H_val] = computeEnergy(obj, alpha)
            alpha_dif = [alpha; obj.alpha_end] - [obj.alpha_start; alpha];
            kappa_eff = alpha_dif ./ obj.eff_lengths - obj.kappa_natural;
            if obj.right_end_free
                kappa_eff(end) = 0;
            end
            E = 0.5 * sum(kappa_eff.^2 .* obj.K_times_l);
            
            if nargout > 1
                tmp = obj.stiffnesses .* kappa_eff;
                grad = tmp(1:end-1) - tmp(2:end);
                
                if nargout > 2
                    E_diag = obj.K_over_l(1:(end-1)) + obj.K_over_l(2:end);
                    if obj.right_end_free
                        E_diag(end) = obj.K_over_l(end-1);
                    end
                    E_offdiag = -obj.K_over_l(2:(end-1));
                    
                    n = obj.num_segments;
                    H_val = [E_diag; E_offdiag; E_offdiag];
                    H_row = [1:n, 1:n-1, 2:n]';
                    H_col = [1:n, 2:n, 1:n-1]';
                end
            end
        end
        
        function [E,grad,H_row,H_col,H_val] = computeGravitationalEnergy(obj, alpha)
            ca = cos(alpha);
            sa = sin(alpha);
            lte = (ca*obj.e_gravity(1) + sa*obj.e_gravity(2)) .* obj.seg_lengths;
            lte_sum = cumsum(lte);
            E = sum(obj.K_times_l(2:end) .* lte_sum);
            
            Kl_upper_sum = flipud(cumsum(flipud(obj.K_times_l(2:end))));
            grad = (-sa*obj.e_gravity(1) + ca*obj.e_gravity(2)) .* obj.seg_lengths .* Kl_upper_sum;
            
            n = obj.num_segments;
            H_row = (1:n)';
            H_col = (1:n)';
            H_val = -lte .* Kl_upper_sum;
        end
        
        % Computes the "constraint energy" sum(lambda_i * C_i); the
        % gradient wrt the primary unknowns; the gradient wrt the Lagrange
        % multipliers (this is just the values of all constraints); the
        % entries into the alpha-alpha part of the Hessian; the entries
        % into the lambda-alpha part of the Hessian (with alpha running
        % along the rows, and lambda running along the columns)
        function [E,grad_alpha,grad_lambda, ...
                Haa_row,Haa_col,Haa_val, ...
                Hal_row,Hal_col,Hal_val] ...
                = computeConstraintEnergy(obj, alpha, lambda)
            
            lca = obj.seg_lengths .* cos(alpha);
            lsa = obj.seg_lengths .* sin(alpha);
            to = [cumsum(lca) cumsum(lsa)];
            grad_lambda = sum(obj.pc_coeff' .* (to(obj.pc_index, :) + obj.p_start'), 2) - obj.pc_value';
            E = dot(grad_lambda, lambda);
            
            if nargout > 1
                nc = length(obj.pc_index);
                grad_ci_alpha = cell(nc, 1);
                n = obj.num_segments;
                grad_alpha = zeros(n, 1);
                for ci=1:nc
                    pi = obj.pc_index(ci);
                    grad_ci_alpha{ci} = obj.pc_coeff(2,ci) * lca(1:pi) - obj.pc_coeff(1,ci) * lsa(1:pi);
                    grad_alpha(1:pi) = grad_alpha(1:pi) + lambda(ci) * grad_ci_alpha{ci};
                end
                
                if nargout > 3
                    pi_max = max(obj.pc_index);
                    pi_sum = sum(obj.pc_index);
                    
                    Haa_val = zeros(pi_max, 1);
                    Haa_row = (1:pi_max)';
                    Haa_col = (1:pi_max)';
                    
                    Hal_val = zeros(pi_sum, 1);
                    Hal_row = zeros(pi_sum, 1);
                    Hal_col = zeros(pi_sum, 1);
                    
                    off = 0;
                    for ci=1:nc
                        pi = obj.pc_index(ci);
                        Haa_val(1:pi) = Haa_val(1:pi) - lambda(ci) * (obj.pc_coeff(1,ci) * lca(1:pi) + obj.pc_coeff(2,ci) * lsa(1:pi));
                        
                        Hal_val(off+1:off+pi) = grad_ci_alpha{ci};
                        Hal_row(off+1:off+pi) = (1:pi)';
                        Hal_col(off+1:off+pi) = ci;
                        
                        off = off + pi;
                    end
                end
            end
        end
        
        function ret = computeAll(obj, alpha, lambda)
            % Compute lagrangian, grad, and hessian for n angles, and
            % Lagr multipliers lambda. The first and last angle are the
            % tangent angles at the first and last vertex. The remaining
            % angles correspond to the linear segments.
            n = obj.num_segments;
            nc = length(obj.pc_index);
            
            [E, Egrad, EH_row, EH_col, EH_val] = obj.computeEnergy(alpha);
            if any(obj.e_gravity ~= 0)
                [U, Ugrad, UH_row, UH_col, UH_val] = obj.computeGravitationalEnergy(alpha);
            else
                U = 0; Ugrad = zeros(n, 1); UH_row = zeros(0,1); UH_col = zeros(0,1); UH_val = zeros(0,1);
            end
            [cE,cgrad_alpha,cgrad_lambda, ...
                Haa_row,Haa_col,Haa_val, ...
                Hal_row,Hal_col,Hal_val] ...
                = obj.computeConstraintEnergy(alpha, lambda);
            
            L = E + cE + U;
            grad =  [cgrad_lambda; Egrad + cgrad_alpha + Ugrad];
            H_val = [EH_val; Haa_val; Hal_val; Hal_val; UH_val];
            H_row = [nc + EH_row; nc + Haa_row; nc + Hal_row; Hal_col; nc + UH_row];
            H_col = [nc + EH_col; nc + Haa_col; Hal_col; nc + Hal_row; nc + UH_col];
            H = sparse(H_row, H_col, H_val, n+nc, n+nc);
            
            ret = struct('L',L,...
                'grad',grad,...
                'H',H,...
                'E',E,...
                'constraints',cgrad_lambda);
        end
        
        function lambda = findOptimalLambda(obj, alpha)
            nc = length(obj.pc_index);
            [~, Egrad] = obj.computeEnergy(alpha);
            [~,~,~, ...
                ~,~,~, ...
                Hal_row,Hal_col,Hal_val] ...
                = obj.computeConstraintEnergy(alpha, zeros(nc,1));
            
            nc = length(obj.pc_index);
            n = obj.num_segments;
            
            V = sparse(Hal_row, Hal_col, Hal_val, n, nc);
            u = Egrad;
            lambda = (V'*V) \ (-(V'*u));
        end
        
        function [alpha, lambda, converged] = optimizeWithNewton(obj, alpha_init, lambda_init, max_iter, grad_norm_threshold)
            num_constr = length(obj.pc_index);
            if nargin < 3
                lambda_init = zeros(num_constr, 1);
            end
            if nargin < 4
                max_iter = 10;
            end
            if nargin < 5
                grad_norm_threshold = 1.e-10;
            end
            
            alpha = alpha_init;
            lambda = lambda_init;
            converged = false;
            for i=1:max_iter
                ret = obj.computeAll(alpha, lambda);
                step = ret.H \ (-ret.grad);
                alpha = alpha + step(num_constr+1:end);
                lambda = lambda + step(1:num_constr);

                grad_norm = norm(ret.grad);
                if grad_norm < grad_norm_threshold
                    converged = true;
                    break;
                end
            end
        end
        
        function Z = solveAdjointJacobiEquation(obj, alpha, lambda, rhs)
            n = obj.num_segments;
            if nargin < 4
                rhs = zeros(n,1);
            end
            ca = cos(alpha);
            sa = sin(alpha);
            Q = sum(-lambda' .* [ca sa], 2);
            Z = zeros(n+2,1);
            for i=n+1:-1:2
                tmp = obj.seg_lengths(i-1)*(rhs(i-1) - Q(i-1)*Z(i)) + obj.K_over_l(i)*(Z(i+1)-Z(i));
                Z(i-1) = Z(i) - tmp / obj.K_over_l(i-1);
            end
        end
        
        function Z = solveJacobiEquation(obj, alpha, lambda, rhs)
            n = obj.num_segments;
            if nargin < 4
                rhs = zeros(n,1);
            end
            ca = cos(alpha);
            sa = sin(alpha);
            Q = sum(-lambda' .* [ca sa], 2);
            Z = zeros(n+2,1);
            Z(2) = obj.eff_lengths(1);
            for i=2:n+1
                tmp = obj.seg_lengths(i-1)*(Q(i-1)*Z(i)-rhs(i-1)) + obj.K_over_l(i-1)*(Z(i)-Z(i-1));
                Z(i+1) = Z(i) + tmp / obj.K_over_l(i);
            end
        end
        
        function [detM, obj_var_lin, obj_var_0] = computeConstrainedStabilityMatrix(obj, alpha, lambda, p)
            xi = obj.solveJacobiEquation(alpha, lambda);
            T1 = -sin(alpha);
            T2 = cos(alpha);
            eta1 = obj.solveJacobiEquation(alpha, lambda, T1);
            eta2 = obj.solveJacobiEquation(alpha, lambda, T2);
            
            xi_bord = xi([1 end]);
            eta1_bord = eta1([1 end]);
            eta2_bord = eta2([1 end]);
            xi = xi(2:end-1);
            eta1 = eta1(2:end-1);
            eta2 = eta2(2:end-1);
            
            M1 = [0;cumsum(T1.*xi.*obj.seg_lengths)];
            M2 = [0;cumsum(T2.*xi.*obj.seg_lengths)];
            N11 = [0;cumsum(T1.*eta1.*obj.seg_lengths)];
            N12 = [0;cumsum(T1.*eta2.*obj.seg_lengths)];
            N21 = [0;cumsum(T2.*eta1.*obj.seg_lengths)];
            N22 = [0;cumsum(T2.*eta2.*obj.seg_lengths)];
            
            M1 = 0.5 * (M1(1:end-1) + M1(2:end));
            M2 = 0.5 * (M2(1:end-1) + M2(2:end));
            N11 = 0.5 * (N11(1:end-1) + N11(2:end));
            N12 = 0.5 * (N12(1:end-1) + N12(2:end));
            N21 = 0.5 * (N21(1:end-1) + N21(2:end));
            N22 = 0.5 * (N22(1:end-1) + N22(2:end));
            
            M = reshape([xi M1 M2 eta1 N11 N21 eta2 N12 N22]', 3, 3, []);
            
            detM = M(1,1,:).*M(2,2,:).*M(3,3,:) + M(1,2,:).*M(2,3,:).*M(3,1,:) + M(1,3,:).*M(2,1,:).*M(3,2,:) ...
                - M(3,1,:).*M(2,2,:).*M(1,3,:) - M(3,2,:).*M(2,3,:).*M(1,1,:) - M(3,3,:).*M(2,1,:).*M(1,2,:);
            detM = detM(:);
            
            if nargout < 2
                return;
            end
            
            % cofactors of M
            pa_xi = p .* reshape(M(2,2,:).*M(3,3,:) - M(2,3,:).*M(3,2,:),[],1);
            pa_eta1 = p .* reshape(M(2,3,:).*M(3,1,:) - M(2,1,:).*M(3,3,:),[],1);
            pa_eta2 = p .* reshape(M(2,1,:).*M(3,2,:) - M(2,2,:).*M(3,1,:),[],1);
            pa_M1 = p .* reshape(M(3,2,:).*M(1,3,:) - M(1,2,:).*M(3,3,:),[],1);
            pa_N11 = p .* reshape(M(1,1,:).*M(3,3,:) - M(1,3,:).*M(3,1,:),[],1);
            pa_N12 = p .* reshape(M(3,1,:).*M(1,2,:) - M(1,1,:).*M(3,2,:),[],1);
            pa_M2 = p .* reshape(M(1,2,:).*M(2,3,:) - M(2,2,:).*M(1,3,:),[],1);
            pa_N21 = p .* reshape(M(2,1,:).*M(1,3,:) - M(1,1,:).*M(2,3,:),[],1);
            pa_N22 = p .* reshape(M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:),[],1);
            
            % Mibar, Nijbar
            M1bar = [0;cumsum(pa_M1 .* obj.seg_lengths)];
            M1bar = M1bar - M1bar(end);
            M1bar = 0.5 * (M1bar(1:end-1) + M1bar(2:end));
            
            M2bar = [0;cumsum(pa_M2 .* obj.seg_lengths)];
            M2bar = M2bar - M2bar(end);
            M2bar = 0.5 * (M2bar(1:end-1) + M2bar(2:end));
            
            N11bar = [0;cumsum(pa_N11 .* obj.seg_lengths)];
            N11bar = N11bar - N11bar(end);
            N11bar = 0.5 * (N11bar(1:end-1) + N11bar(2:end));
            
            N12bar = [0;cumsum(pa_N12 .* obj.seg_lengths)];
            N12bar = N12bar - N12bar(end);
            N12bar = 0.5 * (N12bar(1:end-1) + N12bar(2:end));
            
            N21bar = [0;cumsum(pa_N21 .* obj.seg_lengths)];
            N21bar = N21bar - N21bar(end);
            N21bar = 0.5 * (N21bar(1:end-1) + N21bar(2:end));
            
            N22bar = [0;cumsum(pa_N22 .* obj.seg_lengths)];
            N22bar = N22bar - N22bar(end);
            N22bar = 0.5 * (N22bar(1:end-1) + N22bar(2:end));
            
            % xibar, etaibar
            xi_rhs = M1bar.*T1 + M2bar.*T2 - pa_xi;
            xibar = obj.solveAdjointJacobiEquation(alpha, lambda, xi_rhs);
            xibar_bord = xibar([1 end]);
            xibar = xibar(2:end-1);
            eta1_rhs = T1.*N11bar + T2.*N21bar - pa_eta1;
            eta1bar = obj.solveAdjointJacobiEquation(alpha, lambda, eta1_rhs);
            eta1bar_bord = eta1bar([1 end]);
            eta1bar = eta1bar(2:end-1);
            eta2_rhs = T1.*N12bar + T2.*N22bar - pa_eta2;
            eta2bar = obj.solveAdjointJacobiEquation(alpha, lambda, eta2_rhs);
            eta2bar_bord = eta2bar([1 end]);
            eta2bar = eta2bar(2:end-1);
            
            % Q1 = -T2, Q2 = T1
            Q1 = -T2;
            Q2 = T1;
            % Build system for adjoint equilibrium equation
            S1 = sum(lambda' .* [Q1 Q2], 2);
            S2 = sum(lambda' .* [-T1 -T2], 2);
            lin_term = xibar .* xi .* S2 ...
                + eta1bar .* (eta1.*S2 - Q1) ...
                + eta2bar .* (eta2.*S2 - Q2) ...
                - (M1bar .*Q1 + M2bar .*Q2) .* xi ...
                - (N11bar.*Q1 + N21bar.*Q2) .* eta1 ...
                - (N12bar.*Q1 + N22bar.*Q2) .* eta2;
            constr_const = xi.*xibar + eta1.*eta1bar + eta2.*eta2bar;
            cc1 = sum(Q1.*constr_const.*obj.seg_lengths);
            cc2 = sum(Q2.*constr_const.*obj.seg_lengths);
            
            % minimize int [0.5 * K *(alphabar')^2 + 0.5 * S1 * alphabar^2 + linterm * alphabar]
            % s.t. int T_i * alphabar + cc1 = 0
            grad = [cc1; cc2; obj.seg_lengths .* lin_term];
            Hdiag = obj.K_over_l(1:end-1) + obj.K_over_l(2:end) + obj.seg_lengths .* S1;
            Hoffdiag = -obj.K_over_l(2:end-1);
            Hborder = obj.seg_lengths .* [T1 T2];
            
            n = obj.num_segments;
            Hvals = [Hdiag; Hoffdiag; Hoffdiag; Hborder(:); Hborder(:)];
            Hrows = [3:n+2, 3:n+1, 4:n+2, 3:n+2, 3:n+2, ones(1,n), 2*ones(1,n)]';
            Hcols = [3:n+2, 4:n+2, 3:n+1, ones(1,n), 2*ones(1,n), 3:n+2, 3:n+2]';
            H = sparse(Hrows, Hcols, Hvals, n+2, n+2);
            adjsol = H\(-grad);
            alphabar = adjsol(3:end);
            
            dalpha = ([alpha; obj.alpha_end] - [obj.alpha_start; alpha]) ./ obj.eff_lengths;
            dalphabar = ([alphabar; 0] - [0; alphabar]) ./ obj.eff_lengths;
            dxi = ([xi; xi_bord(2)] - [xi_bord(1); xi]) ./ obj.eff_lengths;
            dxibar = ([xibar; xibar_bord(2)] - [xibar_bord(1); xibar]) ./ obj.eff_lengths;
            deta1 = ([eta1; eta1_bord(2)] - [eta1_bord(1); eta1]) ./ obj.eff_lengths;
            deta1bar = ([eta1bar; eta1bar_bord(2)] - [eta1bar_bord(1); eta1bar]) ./ obj.eff_lengths;
            deta2 = ([eta2; eta2_bord(2)] - [eta2_bord(1); eta2]) ./ obj.eff_lengths;
            deta2bar = ([eta2bar; eta2bar_bord(2)] - [eta2bar_bord(1); eta2bar]) ./ obj.eff_lengths;
            obj_var_lin = dalpha .* dalphabar + dxi .* dxibar + deta1 .* deta1bar + deta2 .* deta2bar;
            obj_var_0 = xibar_bord(1) + eta1bar_bord(1) + eta2bar_bord(1);
        end
        
        function p = computePoints(obj, alpha)
            lca = obj.seg_lengths .* cos(alpha);
            lsa = obj.seg_lengths .* sin(alpha);
            x = [obj.p_start(1), obj.p_start(1) + cumsum(lca')];
            y = [obj.p_start(2), obj.p_start(2) + cumsum(lsa')];
            p = [x;y];
        end
        
        function p = findInflectionPoints(obj, alpha)
            curv = ([alpha; obj.alpha_end] - [obj.alpha_start; alpha]) ./ obj.eff_lengths;
            curv_pos = curv > 0;
            seg_infl = find(curv_pos(1:end-1) ~= curv_pos(2:end));
            sk = (curv(seg_infl) .* obj.seg_lengths(seg_infl)) ./ (curv(seg_infl) - curv(seg_infl+1));
            p = obj.computePoints(alpha);
            p = p(:,seg_infl) + sk' .* [cos(alpha(seg_infl)) sin(alpha(seg_infl))]';
        end
        
        function [line_lin, line_const] = findInflectionPointLine(obj, alpha)
            % Input: Equilibrium curve with at least 2 inflection points
            % Output: <p, line_lin> + line_const = 0 for all inflection
            % points
            curv = ([alpha; obj.alpha_end] - [obj.alpha_start; alpha]) ./ obj.eff_lengths;
            curv_pos = curv > 0;
            seg_infl = find(curv_pos(1:end-1) ~= curv_pos(2:end));
            sk = (curv(seg_infl) .* obj.seg_lengths(seg_infl)) ./ (curv(seg_infl) - curv(seg_infl+1));
            p = obj.computePoints(alpha);
            p = p(:,seg_infl) + sk' .* [cos(alpha(seg_infl)) sin(alpha(seg_infl))]';
            to = p(:,2) - p(:,1);
            dir = to / norm(to);
            line_lin = [-dir(2); dir(1)];
            line_const = -dot(p(:,1), line_lin);
        end
        
        function draw(obj, alpha, color, linespec, linewidth)
            if nargin < 3
                color = [0 0 0];
            end
            if nargin < 4
                linespec = '-';
            end
            if nargin < 5
                linewidth = 0.5;
            end
            p = obj.computePoints(alpha);
            plot(p(1,:), p(2,:), linespec, 'Color', color, 'LineWidth', linewidth);
        end
        
        function exportSVG(obj, filename, scale, min_width, extrude_length)
            hw = obj.stiffnesses' * (0.5 * min_width / min(obj.stiffnesses));
            s = [0;cumsum(obj.seg_lengths)]';
            x = [-extrude_length, s, s(end)+extrude_length];
            x = [x, fliplr(x), x(1)];
            y = [hw(1), hw, hw(end)];
            y = [y, -fliplr(y), y(1)];
            p = [x;y];
            
            % post-process for output
            p(2,:) = -p(2,:);
            p = p - min(p,[],2);
            sz = max(p,[],2);
            
            pstr = sprintf('%.4f,%.4f ',p*scale);
            pstr = pstr(1:end-1);
            
            file = fopen(filename, 'w');
            fprintf(file,'<svg width=\"%.4f\" height=\"%.4f\" xmlns=\"http://www.w3.org/2000/svg\">\n', sz(1)*scale, sz(2)*scale);
            fprintf(file,'<polyline points=\"');
            fprintf(file,pstr);
            fprintf(file,'\" stroke=\"black\" fill=\"none\"/>\n');
            fprintf(file,'</svg>');
            fclose(file);
        end
    end
    
    methods(Static)
        function [alpha, alpha_start, alpha_end] = computeAngles(gamma, a0, a1)
            p_to = gamma(:,2:end) - gamma(:,1:end-1);
            alpha = atan2(p_to(2,:), p_to(1,:))';
            for i=2:length(alpha)
                alpha(i) = fixAngle(alpha(i), alpha(i-1));
            end
            alpha_start = fixAngle(a0, alpha(1));
            alpha_end = fixAngle(a1, alpha(end));
        end
    end
end