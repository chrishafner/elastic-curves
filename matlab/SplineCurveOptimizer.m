classdef SplineCurveOptimizer < handle
    properties
        spline
        
        % Linear constraints
        
        num_lin_constr_eq = 0 % number of linear constraint equations
        lin_constr = cell(1,0) % cell array of strings; supported: t,p (for tangents and positional coordinates of endpoints)
        lin_constr_opt = cell(1,0) % cell array of linear constraint options
        % for tangent constraints: {t, n}, <gamma'(t), n> = 0
        % for positions: {t, e, c}, <gamma(t), e> - c = 0
        lin_constr_grad % constraint gradients (wrt control point coordinates)
        
        % Non-linear constraints
        
        
        num_nonl_constr_eq = 0 % number of non-linear constraint equations
        nonl_constr = cell(1,0) % cell array of strings; supports 'al' (constant arc-lenght), 'infl-line' (all inflections stay on constant line)
        nonl_constr_opt = cell(1,0) % cell array of non-linear constraint options
        % for al: {a}, int ||gamma'(t)|| - a = 0
        % for infl-line: {num_infl, infl_lin, infl_const}, <p, infl_lin> + infl_const = 0 for inflection points p
        
        
        % Inflection points p must satisfy
        % <infl_line_lin, p> + infl_line_const = 0
        % These values are recomputed after all constraints have been
        % satisfied
        infl_line_lin
        infl_line_const
        
        bump_size = 0.1
    end
    
    methods
        function obj = SplineCurveOptimizer(spline)
            obj.spline = spline;
            obj.lin_constr_grad = zeros(0, 2*obj.spline.ncp);
            %obj.addAllStandardLinearConstraints();
            %obj.addInflectionPointLineConstraint();
            %obj.addArcLengthConstraint();
        end
        
        function cp_iter = optimizeStability(obj, stability_goal, r_max, p, max_step_size, tess, max_iter, plot_iterations)
            if nargin < 7
                max_iter = 500;
            end
            if nargin < 8
                plot_iterations = true;
            end
            
            obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
            
            cp_record = Array();
            
            if plot_iterations
                figure(1);
                clf;
                hold on;

                figure(2);
                clf;
                hold on;
            end
            
            p_spline_init = obj.spline.evaluate(linspace(0,obj.spline.t_max,1e3));
            
            success = false;
            for i=1:max_iter
                cp_record.append(obj.spline.cp(:));
                
                [~, dG_dq] = obj.evaluateConstraints();
                dq_dy = null(dG_dq);
                
                dK_dq = obj.computeStiffnessSD(tess)';
                dK_dy = dK_dq * dq_dy;
                
                [elastica, alpha] = obj.spline.convertToElastica(tess);
                
                [r, dr_dK] = SplineCurveOptimizer.softMaxOverMin(elastica.stiffnesses, elastica.eff_lengths, p);
                dr_dK = dr_dK';
                dr_dy = dr_dK * dK_dy;
                
                fprintf("Iteration %u:\n", i);
                fprintf("Ratio: %8f, Soft ratio: %8f\n", max(elastica.stiffnesses) / min(elastica.stiffnesses), r);
                
                lambda = elastica.findOptimalLambda(alpha);

                
                detM = elastica.computeConstrainedStabilityMatrix(alpha, lambda);
                s_mid = elastica.sMid();
                
                % Find conjugate point
                s_conj = SplineCurveOptimizer.findFirstZeroCrossing(s_mid, detM, 0.1);
                if isinf(s_conj)
                    fprintf("No conjugate points; stability at %.8f\n", detM(end));
                    if detM(end) > stability_goal
                        success = true;
                        break;
                    end
                    s_conj = s_mid(end);
                else
                    fprintf("Conjugate point at %.8f / %.8f\n", s_conj, s_mid(end));
                end
                bump = SplineCurveOptimizer.makeBump(s_conj, obj.bump_size, s_mid);
                
                if plot_iterations
                    figure(1);
                    if mod(i,200)==0
                        clf;
                        plot(p_spline_init(1,:), p_spline_init(2,:));
                        hold on;
                    end
                    p_spline = obj.spline.evaluate(linspace(0,obj.spline.t_max,1e3));
                    plot(p_spline(1,:), p_spline(2,:));
                    axis tight equal;
                    
                    figure(2);
                    if mod(i,200)==0
                        clf;
                        hold on;
                    end
                    plot(s_mid, detM);
                end
                
                [~, df_dK_lin, df_dK_bdry0] = elastica.computeConstrainedStabilityMatrix(alpha, lambda, bump);
                df_dK_lin = df_dK_lin';

                df_dy = sum((elastica.eff_lengths' .* df_dK_lin) * dK_dy, 1) + df_dK_bdry0 * dK_dy(1,:);
                dir_y = df_dy' / norm(df_dy);
                step_y = max_step_size * dir_y;
                step_r_1o = dot(dr_dy, step_y);
                
                if r + step_r_1o > r_max
                    % Want <dr_dy, step_y> + (r-r_max) = 0
                    % Divide by ||dr_dq|| to get <a, step_y> + b = 0
                    dr_dy_len = norm(dr_dy);
                    a = dr_dy' / dr_dy_len;
                    b = (r - r_max) / dr_dy_len;
                    step_y = SplineCurveOptimizer.clampStep(a,b,step_y);
                end
                
                step_r_1o = dot(dr_dy, step_y);
                fprintf('Predicted ratio change: %8f\n', step_r_1o);
                
                step_q = dq_dy * step_y;
                
                cp_new = obj.spline.cp + reshape(step_q,2,[]);
                obj.spline.cp = cp_new;
                
                ok = obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
                if ~ok
                    break;
                end
            end
            
            if success
                fprintf("Success after %u iterations!\n", i);
            else
                fprintf("Failure after %u iterations!\n", i);
            end
            
            % file_out = fopen("spline_output.txt","w");
            % fprintf(file_out, "%u %u\n", obj.spline.degree, obj.spline.ncp);
            % fprintf(file_out, "%.15e %.15e\n", obj.spline.cp);
            % fclose(file_out);
            
            cp_iter = reshape(cp_record.get(), 2, obj.spline.ncp, []);
        end
        
        function cp_iter = optimizeStiffness(obj, target_ratio, p, max_step_size, tess)
            if nargin < 4
                max_step_size = 0.1;
            end
            if nargin < 5
                tess = 20;
            end
            
            obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
            
            max_iter = 5000;
            
            cp_record = Array();
            cp0 = obj.spline.cp;
            for i=1:max_iter+1
                cp_record.append(obj.spline.cp(:));
                
                % Compute orthogonal subspace of constraint tangents
                [~, dG_dq] = obj.evaluateConstraints(); %dG/dq
                dq_dy = null(dG_dq); % dq/dy
                
                % Compute objective and gradient
                [K, ~, weights] = obj.computeStiffness(tess);
                dK_dq = obj.computeStiffnessSD(tess)'; % dK/dq
                dK_dy = dK_dq * dq_dy; % dK/dy
                [f, df_dK] = obj.softMaxOverMin(K,weights,p); %f, df/dK
                df_dy = df_dK * dK_dy; % df/dy
                df_dy_len = norm(df_dy);
                
                fprintf("Iteration %u\n", i-1);
                fprintf("Soft ratio: %8f, real ratio: %8f\n", f, max(K)/min(K));
                if i==max_iter+1 || f <= target_ratio+1.e-4
                    break;
                end
                
                % Compute search direction in y-space
                dir_y = df_dy' / df_dy_len;
                % Compute step size that reaches goal to first order
                h_first_order = (target_ratio - f) / df_dy_len;
                h = h_first_order;
                % Clamp step size
                if abs(h) > max_step_size
                    h = max_step_size * sign(h_first_order);
                end
                
                % Take step in q-space and take it
                step_y = h * dir_y;
                step_q = dq_dy * step_y;
                obj.spline.cp = obj.spline.cp + reshape(step_q, 2, []);
                obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
            end
            cp_dist = sqrt(sum((obj.spline.cp - cp0).^2, [1 2]));
            fprintf("CP dist: %.4f\n", cp_dist);
            
            cp_iter = reshape(cp_record.get(), 2, obj.spline.ncp, []);
        end
        
        function cp_iter = optimizeShape(obj, cp_goal, r_max, p, sigma_min, rate, tess, max_iter)
            if nargin < 7
                max_iter = 300;
            end
            obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
            
            cp_record = Array();
            for i=1:max_iter+1
                cp_record.append(obj.spline.cp(:));
                 % Compute orthogonal subspace of constraint tangents
                [~, dG_dq] = obj.evaluateConstraints(); %dG/dq
                dq_dy = null(dG_dq); % dq/dy
                
                % Compute stiffness gradient
                dK_dq = obj.computeStiffnessSD(tess)'; % dK/dq
                dK_dy = dK_dq * dq_dy; % dK/dy
                
                % Compute ratio gradient
                [elastica, alpha] = obj.spline.convertToElastica(tess);
                [r, dr_dK] = SplineCurveOptimizer.softMaxOverMin(elastica.stiffnesses, elastica.eff_lengths, p);
                dr_dK = dr_dK';
                dr_dy = dr_dK * dK_dy;
                
                
                lambda = elastica.findOptimalLambda(alpha);
                s_mid = elastica.sMid();
                bump = SplineCurveOptimizer.makeBump(s_mid(end), 0.1, s_mid);
                bump_vol = sum(elastica.seg_lengths .* bump);
                [detM, dsigma_dK_lin, dsigma_dK_bdry0] = elastica.computeConstrainedStabilityMatrix(alpha, lambda, bump);
                dsigma_dK_lin = dsigma_dK_lin';
                
                sigma = sum(elastica.seg_lengths .* bump .* detM) / bump_vol;
                dsigma_dy = (sum((elastica.eff_lengths' .* dsigma_dK_lin) * dK_dy, 1) + dsigma_dK_bdry0 * dK_dy(1,:)) / bump_vol;
                
                fprintf("Iteration %u\n", i-1);
                fprintf("Soft ratio: %8f, real ratio: %8f\n", r, max(elastica.stiffnesses)/min(elastica.stiffnesses));
                fprintf("Sigma: %8f\n", sigma);
                
                df_dq = -reshape(obj.spline.cp - cp_goal,[],1);
                df_dy = dq_dy' * df_dq;
                step_y = rate * df_dy;
                step_r_1o = dot(dr_dy, step_y);
                step_sigma_1o = dot(dsigma_dy, step_y);
                if r + step_r_1o > r_max || sigma + step_sigma_1o < sigma_min
                    % Want <dr_dy, step_y> + (r-r_max) = 0
                    % Divide by ||dr_dq|| to get <a, step_y> + b = 0
                    dr_dy_len = norm(dr_dy);
                    a1 = -dr_dy' / dr_dy_len;
                    b1 = (r_max - r) / dr_dy_len;
                    dsigma_dy_len = norm(dsigma_dy);
                    a2 = dsigma_dy' / dsigma_dy_len;
                    b2 = (sigma - sigma_min) / dsigma_dy_len;
                    % step_y = SplineCurveOptimizer.clampStep(a1,b1,step_y);
                    step_y = SplineCurveOptimizer.projectOntoIntersectionOfHalfSpaces(a1,b1,a2,b2,step_y);
                end
                
                step_r_1o = dot(dr_dy, step_y);
                fprintf('Obj: %.8f\n', sum(sum((obj.spline.cp - cp_goal).^2)));
%                 fprintf('Predicted ratio change: %8f\n', step_r_1o);
                
                step_q = dq_dy * step_y;
                
                cp_new = obj.spline.cp + reshape(step_q,2,[]);
                obj.spline.cp = cp_new;
                ok = obj.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
                if ~ok
                    break;
                end
            end
            cp_iter = reshape(cp_record.get(), 2, obj.spline.ncp, []);
        end
        
        function [ret, t, weights] = computeStiffness(obj, tess)
            if nargin < 2
                tess = 20;
            end
            t = linspace(0, obj.spline.t_max, tess*obj.spline.t_max+1);
            p = obj.spline.evaluate(t);
            curv = obj.spline.curvature(t);
            ret = (sum(p .* obj.infl_line_lin, 1) + obj.infl_line_const) ./ curv;
            if nargout > 2
                seg_lengths = sqrt(sum((p(:,2:end) - p(:,1:end-1)).^2, 1));
                weights = 0.5*[seg_lengths(1), (seg_lengths(1:end-1) + seg_lengths(2:end)), seg_lengths(end)];
            end
        end
        
        function ret = computeStiffnessSD(obj, tess)
            if nargin < 2
                tess = 20;
            end
            t = linspace(0, obj.spline.t_max, tess*obj.spline.t_max+1);
            
            cpi_first = obj.spline.firstControlPointIndexAt(t);
            p = obj.spline.evaluate(t);
            curv = obj.spline.curvature(t);
            num = sum(p .* obj.infl_line_lin, 1) + obj.infl_line_const;
            
            curv_sd = obj.spline.curvatureSD(t);
            w_sd = obj.spline.evaluateWeights(t);
            num_sd = zeros(size(curv_sd));
            num_sd(1:2:end,:) = w_sd * obj.infl_line_lin(1);
            num_sd(2:2:end,:) = w_sd * obj.infl_line_lin(2);
            
            grad = num_sd ./ curv - (num .* curv_sd) ./ (curv.^2);
            ret = zeros(2*obj.spline.ncp, length(t));
            blocks = [1, find(cpi_first(1:end-1) ~= cpi_first(2:end)) + 1, length(t)+1];
            nr = 2*(obj.spline.degree+1);
            for i=1:length(blocks)-1
                ret(2*i-1:2*i-2+nr, blocks(i):blocks(i+1)-1) = grad(:,blocks(i):blocks(i+1)-1);
            end
        end
        
        function ok = satisfyConstraintsWithUnderdeterminedNewton(obj, threshold, max_steps)
            if nargin < 3
                max_steps = 10;
            end
            ok = false;
            for i=1:max_steps
                [f, grad_f] = obj.evaluateConstraints();
                if any(isinf(f))
                    break;
                end
                if norm(f) < threshold
                    ok = true;
                    break;
                end
                x = obj.spline.cp(:);
                step = lsqminnorm(grad_f, -f);
                x_new = x + step;
                obj.spline.cp = reshape(x_new, 2, []);
            end
            [obj.infl_line_lin, obj.infl_line_const] = obj.computeInflectionLine();
        end
        
        %%%%%%%%%%%%%%%% Linear Constraints %%%%%%%%%%%%%%%%%%%%%
        
        function addAllStandardLinearConstraints(obj)
            % Add constraints to fix tangent and position at the two
            % endpoints
            for t = [0 obj.spline.t_max]
                obj.addTangentConstraint(t);
                obj.addPositionConstraint(t, [1; 0]);
                obj.addPositionConstraint(t, [0; 1]);
            end
        end
        
        function addTangentConstraint(obj, t, tangent)
            % Fix direction of gamma'(t)'
            obj.lin_constr = [obj.lin_constr, 't'];
            
            if nargin < 3
                tangent = obj.spline.evaluateD(t);
            end
            normal = [-tangent(2); tangent(1)];
            normal = normal / norm(normal);
            obj.lin_constr_opt = [obj.lin_constr_opt, {{t, normal}}];
            
            dw = obj.spline.evaluateWeightsD(t);
            cpi_all = obj.spline.firstControlPointIndexAt(t) + (0:obj.spline.degree)';
            grad = zeros(1, 2*obj.spline.ncp);
            grad(cpi_all*2-1) = dw * normal(1);
            grad(cpi_all*2) = dw * normal(2);
            obj.lin_constr_grad = [obj.lin_constr_grad; grad];
            
            obj.num_lin_constr_eq = obj.num_lin_constr_eq + 1;
        end
        
        function addPositionConstraint(obj, t, e)
            % Fix <gamma(t), e>
            obj.lin_constr = [obj.lin_constr, 'p'];
            
            pos = obj.spline.evaluate(t);
            c = dot(pos, e);
            obj.lin_constr_opt = [obj.lin_constr_opt, {{t, e, c}}];
            
            w = obj.spline.evaluateWeights(t);
            cpi_all = obj.spline.firstControlPointIndexAt(t) + (0:obj.spline.degree)';
            grad = zeros(1, 2*obj.spline.ncp);
            grad(cpi_all*2-1) = w * e(1);
            grad(cpi_all*2) = w * e(2);
            obj.lin_constr_grad = [obj.lin_constr_grad; grad];
            
            obj.num_lin_constr_eq = obj.num_lin_constr_eq + 1;
        end
        
        function ret = evaluateLinearConstraints(obj)
            nc = length(obj.lin_constr);
            ret = zeros(obj.num_lin_constr_eq, 1);
            off = 0;
            for ci=1:nc
                constr_type = obj.lin_constr{ci};
                switch constr_type
                    case 't'
                        t = obj.lin_constr_opt{ci}{1};
                        normal_0 = obj.lin_constr_opt{ci}{2};
                        tangent = obj.spline.evaluateD(t);
                        val = dot(tangent, normal_0);
                        ret(off+1) = val;
                        off = off + 1;
                    case 'p'
                        t = obj.lin_constr_opt{ci}{1};
                        e = obj.lin_constr_opt{ci}{2};
                        c = obj.lin_constr_opt{ci}{3};
                        p = obj.spline.evaluate(t);
                        val = dot(p,e) - c;
                        ret(off+1) = val;
                        off = off + 1;
                end
            end
        end
        
        %%%%%%%%%%%%%%%% Non-linear Constraints %%%%%%%%%%%%%%%%%%%%%
        
        function addArcLengthConstraint(obj)
            obj.nonl_constr = [obj.nonl_constr, 'al'];
            
            al = obj.spline.arcLength();
            obj.nonl_constr_opt = [obj.nonl_constr_opt, {{al}}];
            
            obj.num_nonl_constr_eq = obj.num_nonl_constr_eq + 1;
        end
        
        function addInflectionPointLineConstraint(obj)
            obj.nonl_constr = [obj.nonl_constr, 'infl-line'];
            
            % Line connecting first and last inflection point
            [~, p_infl] = obj.spline.findInflectionPoints();
            num_infl = size(p_infl, 2);
            [infl_lin, infl_const] = obj.computeInflectionLine();
            
            obj.nonl_constr_opt = [obj.nonl_constr_opt, {{num_infl, infl_lin, infl_const}}];
            
            obj.num_nonl_constr_eq = obj.num_nonl_constr_eq + num_infl;
        end

%         function addInflectionPointCollinearConstraint(obj)
%             obj.nonl_constr = [obj.nonl_constr, 'infl-collin'];
%             
%             t_infl = obj.spline.findInflectionPoints();
%             num_infl = length(t_infl);
%             obj.nonl_constr_opt = [obj.nonl_constr_opt, {{num_infl}}];
%             
%             obj.num_nonl_constr_eq = obj.num_nonl_constr_eq + 1;
%         end
        
        function [infl_lin, infl_const] = computeInflectionLine(obj)
            % Line connecting first and last inflection point
            [~, p_infl] = obj.spline.findInflectionPoints();
            to = p_infl(:,end) - p_infl(:,1);
            dir = to / norm(to);
            infl_lin = [-dir(2); dir(1)];
            infl_const = -dot(infl_lin, p_infl(:,1));
            
            % Check if sign of line equation needs to be flipped
            t_samp = linspace(0, obj.spline.t_max, 50);
            curv_samp = obj.spline.curvature(t_samp);
            [~,i_max] = max(abs(curv_samp));
            p = obj.spline.evaluate(t_samp(i_max));
            ret = (dot(p, infl_lin) + infl_const) / curv_samp(i_max);
            if ret < 0
                infl_lin = -infl_lin;
                infl_const = -infl_const;
            end
        end
        
        function [ret, grad]  = evaluateNonlinearConstraints(obj)
            nc = length(obj.nonl_constr);
            neq = obj.num_nonl_constr_eq;
            
            ret = zeros(neq, 1);
            grad = zeros(neq, 2 * obj.spline.ncp);
            off = 0;
            
            for ci=1:nc
                constr_type = obj.nonl_constr{ci};
                switch constr_type
                    case 'al'
                        al = obj.spline.arcLength();
                        al0 = obj.nonl_constr_opt{ci}{1};
                        ret(off+1) = al - al0;
                        if nargout > 1
                            dal = obj.spline.arcLengthSD();
                            grad(off+1,:) = dal(:)';
                        end
                        off = off + 1;
                    case 'infl-line'
                        num_infl = obj.nonl_constr_opt{ci}{1};
                        infl_lin = obj.nonl_constr_opt{ci}{2};
                        infl_const = obj.nonl_constr_opt{ci}{3};
                        [t_infl, p_infl] = obj.spline.findInflectionPoints();
                        if length(t_infl) ~= num_infl
                            ret(off+(1:num_infl)) = inf;
                            grad(off+(1:num_infl), :) = nan;
                        else
                            ret(off+(1:num_infl)) = sum(p_infl .* infl_lin, 1)' + infl_const;
                            if nargout > 1
                                cpi = obj.spline.firstControlPointIndexAt(t_infl);
                                dofi = cpi*2-1;
                                for i=1:num_infl
                                    [~, dp_dq] = obj.spline.inflectionPointSD(t_infl(i));
                                    dconstr_dq = sum(dp_dq.*infl_lin,1);
                                    grad(off+i, dofi(i):dofi(i)+(2*obj.spline.degree+1)) = dconstr_dq;
                                end
                            end
                        end
                        off = off + num_infl;
%                     case 'infl-collin'
%                         num_infl = obj.nonl_constr_opt{ci}{1};
%                         [t_infl, p_infl] = obj.spline.findInflectionPoints();
%                         if length(t_infl) ~= num_infl
%                             ret(off+1) = inf;
%                             grad(off+1, :) = nan;
%                         else
%                             to = p_infl(:,end) - p_infl(:,1);
%                             dir = to / norm(to);
%                             infl_lin = [-dir(2); dir(1)];
%                             infl_const = -dot(infl_lin, p_infl(:,1));
%                         end
                end
            end
        end
        
        function [ret, grad] = evaluateConstraints(obj)
            if nargout > 1
                [ret_nonl, grad_nonl] = obj.evaluateNonlinearConstraints();
                ret = [obj.evaluateLinearConstraints(); ret_nonl];
                grad = [obj.lin_constr_grad; grad_nonl];
            else
                ret_nonl = obj.evaluateNonlinearConstraints();
                ret = [obj.evaluateLinearConstraints(); ret_nonl];
            end
        end
    end
    
    methods(Static)
        function bump = makeBump(center, half_width, s)
            x = (s-center)/half_width;
            bump = zeros(size(x));
            bump(abs(x)<1) = exp(-1./(1-x(abs(x)<1).^2));
        end
        
        function [obj, grad] = softMax(K,w,p)
            k = max(K);
            KK = K / k;
            obj = k * sum(KK.^p .* w) .^ (1/p);
            if nargout > 1
                grad = (obj / k)^(1-p) * (KK.^(p-1) .* w);
            end
        end
        
        function [obj, grad] = softMaxOverMin(K,w,p)
            if nargout > 1
                [mx, mx_grad] = SplineCurveOptimizer.softMax(K,w,p);
                [mx_inv, mx_inv_grad] = SplineCurveOptimizer.softMax(1./K,w,p);
                mx_inv_grad = -mx_inv_grad ./ (K.^2);
                obj = mx * mx_inv;
                grad = mx_grad * mx_inv + mx * mx_inv_grad;
            else
                mx = SplineCurveOptimizer.softMax(K,w,p);
                mx_inv = SplineCurveOptimizer.softMax(1./K,w,p);
                obj = mx * mx_inv;
            end
        end
        
        function x_zero = findFirstZeroCrossing(x, f, frac_offset)
            off = floor(length(x) * frac_offset);
            i = find(sign(f(off:end-1)) ~= sign(f(off+1:end)), 1) + off;
            if isempty(i)
                x_zero = inf;
            else
                t = f(i-1) / (f(i-1) - f(i));
                x_zero = x(i-1)*(1-t) + x(i)*t;
            end
        end
        
        % Finds an x closest to x0 such that <x,a>+b=0 and ||x||=||x0||
        % If it does not exist, finds x closest to {x:<x,a>+b=0} and such
        % that ||x||=||x0||
        
        % Input must satisfy ||a||=1
        function x = clampStep(a,b,x0)
            z = x0 - dot(x0,a)*a;
            z_sq_len = dot(z,z);
            
            while z_sq_len < 1.e-20
                z = randn(size(x0));
                z = z - dot(z,a)*a;
                z_sq_len = dot(z,z);
            end
                
            c_sq = (dot(x0,x0) - b*b) / z_sq_len;
            if c_sq < 0
                x = -sign(b) * norm(x0) * a;
                return;
            end
            c_pm = sqrt(c_sq);
            x1 = -b*a + c_pm*z;
            x2 = -b*a - c_pm*z;
            x1_sq_dist = sum((x1 - x0).^2);
            x2_sq_dist = sum((x2 - x0).^2);
            if x1_sq_dist < x2_sq_dist
                x = x1;
            else
                x = x2;
            end
        end
        
        % Finds x which is closest to x0 and satisfies <ai,x>+bi >= 0
        % Input must satisfy ||ai||=1
        function x = projectOntoIntersectionOfHalfSpaces(a1,b1,a2,b2,x0)
            r1 = -b1-dot(x0,a1);
            r2 = -b2-dot(x0,a2);
            % Projection onto the two planes
            x_p1 = x0 + r1*a1;
            x_p2 = x0 + r2*a2;
            a = dot(a1,a2);
            if abs(1-a) < 1.e-10
                x = x_p1;
                return;
            end
            
            % Projection onto the intersection line
            fac = 1/(1-a*a);
            x_pl = x0 + (fac*(r1-a*r2))*a1 + (fac*(r2-a*r1))*a2;
            
            % Sort the three candidates by distance to x0, and return the
            % closest feasible point
            x_cand = [x_p1, x_p2, x_pl];
            dist_sq = sum((x_cand - x0).^2, 1);
            ok = [dot(x_p1,a2)+b2 >= 0, dot(x_p2,a1)+b1 >= 0, true];
            [~,I] = sort(dist_sq);
            if ok(I(1))
                x = x_cand(:,I(1));
            elseif ok(I(2))
                x = x_cand(:,I(2));
            else
                x = x_cand(:,I(3));
            end
        end
    end
end