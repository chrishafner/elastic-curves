classdef SplineCurveTools < handle
    methods(Static)
        function cp_final = enforceInflectionConstraint(spline, line_lin, line_const)
            % Enforces for all inflection points p:
            % <p,line_lin> + line_const = 0.
            
            options = optimoptions('fmincon', ...
                'Algorithm', 'sqp', ...
                'ConstraintTolerance', 1.e-10, ...
                'Display', 'off', ... %'iter-detailed', ...
                'SpecifyConstraintGradient', true, ...
                'SpecifyObjectiveGradient', true);
            
            t_infl = spline.findInflectionPoints();
            num_infl = length(t_infl);
            
            x0 = spline.cp(:);
            obj_fun = @(x) SplineCurveTools.controlPointObjective(x0,x);
            nonlcon_fun = @(x) SplineCurveTools.infl_constr(spline, num_infl, x, line_lin, line_const);
            
            [x_final,fval,exitflag,output] = fmincon(obj_fun, x0, [], [], [], [], [], [], nonlcon_fun, options);
            
            cp_final = reshape(x_final, 2, []);
        end
        
        % endpoints_fixed is a boolean 4-vector [p0x,p0y,p1x,p1y]
        function [ret, iterations, iterations_obj_val] = optimizeStiffnessRatio(spline, R_target, t_tol_max, t_tol_min, endpoints_fixed, plot_iterations, p_approx, target_change, auto_remove_inflections, tangents_fixed)
            if nargin < 3
                t_tol_max = 0.25;
                t_tol_min = 0.05;
            end
            if nargin < 5
                endpoints_fixed = [false;false;false;false];
            end
            if nargin < 6
                plot_iterations = false;
            end
            if nargin < 7
                p_approx = 40;
            end
            if nargin < 8
                target_change = 0.05;
            end
            if nargin < 9
                auto_remove_inflections = true;
            end
            if nargin < 10
                tangents_fixed = [false;false];
            end
            
            
            % Build linear endpoint constraints
            [lin_constr_mat, lin_constr_rhs] = SplineCurveTools.buildEndpointConstraints(spline, endpoints_fixed, tangents_fixed);
            
            % Build preconditioning matrix for step direction
            t_dense = linspace(0,spline.t_max,2e2);
            cp_initial = spline.cp;
            gamma_initial = spline.evaluate(t_dense);
            BtB = SplineCurveTools.buildPreconditioningMatrix(spline);
            
            if auto_remove_inflections
                [t_infl, p_infl] = spline.findInflectionPoints([t_tol_min, t_tol_min]);
                t_threshold = 0.05 * spline.t_max;
                if length(t_infl) > 2
                    [~, I] = sort(abs(t_infl - spline.t_max*0.5));
                    t_keep = t_infl(I([1 2]));
                    t_remove = t_infl(I([3 4]));
                    t_remove_dist = max(min([t_remove'; spline.t_max - t_remove'], [], 1));
                    t_keep_dist = min(min([t_keep'; spline.t_max - t_keep'], [], 1));
                    t_threshold = 0.5 * (t_remove_dist + t_keep_dist);
                end
                SplineCurveTools.removeInflectionPoints(spline, endpoints_fixed, t_threshold, t_tol_max, plot_iterations);
            end
            
            cp_preprocessed = spline.cp;
            
            if plot_iterations
                figure(1);
                clf;
                hold on;

                figure(2);
                clf;
                hold off;
            end
            
            tol_infl = [true true];
            [t_infl, ~] = spline.findInflectionPoints([t_tol_min, t_tol_min]);
            if any(t_infl < spline.t_max*0.5)
                tol_infl(1) = true;
            end
            if any(t_infl > spline.t_max*0.5)
                tol_infl(2) = true;
            end
            
            iterations = CellArray();
            iterations_obj_val = Array();
            
            K_start = [];
            max_iter = 1e3;
            for i=1:max_iter
                gamma = spline.evaluate(t_dense);
                
                % Look for inflection points
                t_tol = [t_tol_min t_tol_min];
                if t_tol_max > t_tol_min
                    if tol_infl(1)
                        t_tol(1) = t_tol_max;
                    end
                    if tol_infl(2)
                        t_tol(2) = t_tol_max;
                    end
                end
                kappa = spline.curvature(t_dense);
                [t_infl, p_infl] = spline.findInflectionPoints(t_tol);
                
                if plot_iterations
                    figure(1);
                    clf;
                    hold on;
                    plot(cp_initial(1,:), cp_initial(2,:), 'Color',[1 0.5 0.5],'LineWidth',3);
                    plot(spline.cp(1,:), spline.cp(2,:), 'Color',[1 0 0]);

                    plot(gamma_initial(1,:), gamma_initial(2,:), 'Color',[0.7 0.7 0.7], 'LineWidth',3);
                    plot(gamma(1,:), gamma(2,:), 'Color',[0 0 0]);

                    scatter(p_infl(1,:), p_infl(2,:), 80, [0 0 1], '.');
                    axis tight equal;
                end
                
                % Cases for different numbers of inflection points
                if length(t_infl) >= 2
                    if plot_iterations
                        figure(1);
                        plot(p_infl(1,:), p_infl(2,:), 'Color',[0 0 1], 'LineWidth',2);
                    end
                    
                    infl_to = p_infl(:,2) - p_infl(:,1);
                    a = [-infl_to(2); infl_to(1)];
                    b = -dot(p_infl(:,1), a);
                    K = (sum(a.*gamma,1)+b) ./ kappa;
                    if K(1) < 0
                        K = -K;
                        a = -a;
                        b = -b;
                    end
                    if any(K <= 0)
                        fprintf('Inflection line crosses curve.');
                    end
                else
                    lpopt = LPStiffnessOptimizer(gamma, kappa);
                    if length(t_infl) == 1
                        lpopt.gamma_infl = p_infl;
                        [~, b, a] = lpopt.optimizeWithInflections();
                    else
                        [~, b, a] = lpopt.optimizeSimple();
                    end
                end
                
                [objval, K, objgrad, ~] = SplineCurveTools.computeSoftR(spline, p_approx, a, b, t_dense, t_tol);
                R_current = max(K/min(K));
                if i==1
                    K_start = K;
                    K_plot_y_max = min(max(K_start/min(K_start)), 3*R_target);
                end
                
                if plot_iterations
                    figure(2);
                    clf;
                    plot(t_dense, K_start/min(K_start), 'LineWidth',3, 'Color',[0.5 0.5 0.5]);
                    hold on;
                    plot(t_dense, K/min(K), 'Color',[0 0 0]);
                    ylim([0,K_plot_y_max]);
                end
                
                fprintf('Obj val: %5.2f;   R: %5.2f\n', objval, R_current);
                
                iterations.append(spline.cp);
                iterations_obj_val.append(objval);
                
                if objval < R_target || R_current < R_target
                    break;
                end
                
                target_change_f = -objval * target_change;
                options = optimoptions('quadprog','Display','none');
                dq = quadprog(BtB, zeros(2*spline.ncp,1),[],[],...
                    [objgrad'; lin_constr_mat],[target_change_f;zeros(size(lin_constr_rhs))],...
                    [],[],[],options);
%                 dq = objgrad / dot(objgrad, objgrad) * target_change_f;
                spline.cp = spline.cp + reshape(dq, 2, []);
                
                tol_infl = [false false];
                if any(t_infl < spline.t_max*0.5)
                    tol_infl(1) = true;
                end
                if any(t_infl > spline.t_max*0.5)
                    tol_infl(2) = true;
                end
            end
            
            minK = min(K);
            a = a / minK;
            b = b / minK;
            K = K / minK;
            ret = struct('cp_initial', cp_initial, 'cp_final', spline.cp, 'cp_preprocessed', cp_preprocessed, 'degree', spline.degree, ...
                'R_soft', objval, 'R_hard', max(K), 'K', K, 'a', b, 'b', a, ...
                'gamma_sampled',spline.evaluate(t_dense),'t_sampled',t_dense);
            iterations = iterations.get();
            iterations_obj_val = iterations_obj_val.get();
        end
        
        % Removes inflection points in (0,dt) and (t_max-dt, t_max) by
        % pushing them out of (-t_tol, t_max+t_tol)
        function removeInflectionPoints(spline, endpoints_fixed, dt, t_tol, plot_iterations)
            if nargin < 5
                plot_iterations = false;
            end
            [lin_constr_mat, lin_constr_rhs] = SplineCurveTools.buildEndpointConstraints(spline, endpoints_fixed);
            BtB = SplineCurveTools.buildPreconditioningMatrix(spline);
            
            if plot_iterations
                figure(1);
                clf;
                hold on;
            
                cp_initial = spline.cp;
                t_dense = linspace(0,spline.t_max,2e2);
                gamma_initial = spline.evaluate(t_dense);
            end
            
            iter_max = 100;
            for i=1:iter_max
                [t_infl, p_infl] = spline.findInflectionPoints([t_tol, t_tol]);
                bool_prob = (-t_tol <= t_infl & t_infl <= dt) | (spline.t_max-dt <= t_infl & t_infl <= spline.t_max+t_tol);
                t_prob = t_infl(bool_prob);
                
                if plot_iterations
                    figure(1);
                    clf;
                    hold on;
                    plot(cp_initial(1,:), cp_initial(2,:), 'Color',[1 0.5 0.5],'LineWidth',3);
                    plot(spline.cp(1,:), spline.cp(2,:), 'Color',[1 0 0]);

                    gamma = spline.evaluate(t_dense);
                    plot(gamma_initial(1,:), gamma_initial(2,:), 'Color',[0.7 0.7 0.7], 'LineWidth',3);
                    plot(gamma(1,:), gamma(2,:), 'Color',[0 0 0]);

                    scatter(p_infl(1,:), p_infl(2,:), 80, [0 0 1], '.');
                    axis tight equal;
                end
                
                if isempty(t_prob)
                    break;
                end
                
                t_prob = t_prob(1);
%                 p_prob = p_prob(:,1);
                move_dir = 1;
                if -t_tol <= t_prob && t_prob <= dt
                    move_dir = -1;
                end
                
                cpi = spline.firstControlPointIndexAt(t_prob);
                dofi = (2*cpi-1):(2*cpi+2*spline.degree);
                dt_infl_dq = spline.inflectionPointSD(t_prob);
                dt_dq = zeros(2*spline.ncp, 1);
                dt_dq(dofi) = dt_infl_dq;
                
                target_change = move_dir * 0.1 * (dt + t_tol);
                options = optimoptions('quadprog','Display','none');
                step = quadprog(BtB, zeros(2*spline.ncp,1),[],[],[dt_dq'; lin_constr_mat],[target_change; zeros(size(lin_constr_rhs))],[],[],[],options);
            
                spline.cp = spline.cp + reshape(step,2,[]);
            end
        end
        
        function removeBoundaryInflectionPoint(spline, endpoints_fixed, remove_last, t_tol, plot_iterations)
            if nargin < 5
                plot_iterations = false;
            end
            [lin_constr_mat, lin_constr_rhs] = SplineCurveTools.buildEndpointConstraints(spline, endpoints_fixed);
            BtB = SplineCurveTools.buildPreconditioningMatrix(spline);
            
            if plot_iterations
                figure(1);
                clf;
                hold on;
            
                cp_initial = spline.cp;
                t_dense = linspace(0,spline.t_max,2e2);
                gamma_initial = spline.evaluate(t_dense);
            end
            
            t_infl = spline.findInflectionPoints([t_tol, t_tol]);
            num_infl_initial = length(t_infl);
            
            iter_max = 100;
            for i=1:iter_max
                [t_infl, p_infl] = spline.findInflectionPoints([t_tol, t_tol]);
                if plot_iterations
                    figure(1);
                    clf;
                    hold on;
                    plot(cp_initial(1,:), cp_initial(2,:), 'Color',[1 0.5 0.5],'LineWidth',3);
                    plot(spline.cp(1,:), spline.cp(2,:), 'Color',[1 0 0]);

                    gamma = spline.evaluate(t_dense);
                    plot(gamma_initial(1,:), gamma_initial(2,:), 'Color',[0.7 0.7 0.7], 'LineWidth',3);
                    plot(gamma(1,:), gamma(2,:), 'Color',[0 0 0]);

                    scatter(p_infl(1,:), p_infl(2,:), 80, [0 0 1], '.');
                    axis tight equal;
                end
                
                num_infl = length(t_infl);
                if num_infl < num_infl_initial
                    break;
                end
                
                t_pick = t_infl(1);
                move_dir = -1;
                if remove_last
                    t_pick = t_infl(end);
                    move_dir = 1;
                end
                
                cpi = spline.firstControlPointIndexAt(t_pick);
                dofi = (2*cpi-1):(2*cpi+2*spline.degree);
                dt_infl_dq = spline.inflectionPointSD(t_pick);
                dt_dq = zeros(2*spline.ncp, 1);
                dt_dq(dofi) = dt_infl_dq(:);
                
                target_change = move_dir * 0.1;
                options = optimoptions('quadprog','Display','none');
                step = quadprog(BtB, zeros(2*spline.ncp,1),[],[],[dt_dq'; lin_constr_mat],[target_change; zeros(size(lin_constr_rhs))],[],[],[],options);
            
                spline.cp = spline.cp + reshape(step,2,[]);
            end
        end
        
        % Cancels the inflection points with indices ind_infl and
        % ind_infl+1
        function removeInflectionPointPair(spline, endpoints_fixed, ind_infl, plot_iterations)
            if nargin < 4
                plot_iterations = false;
            end
            [lin_constr_mat, lin_constr_rhs] = SplineCurveTools.buildEndpointConstraints(spline, endpoints_fixed);
            BtB = SplineCurveTools.buildPreconditioningMatrix(spline);
            
            if plot_iterations
                figure(1);
                clf;
                hold on;
            
                cp_initial = spline.cp;
                t_dense = linspace(0,spline.t_max,2e2);
                gamma_initial = spline.evaluate(t_dense);
            end
            
            t_infl = spline.findInflectionPoints();
            num_infl_initial = length(t_infl);
            
            iter_max = 100;
            for i=1:iter_max
                [t_infl, p_infl] = spline.findInflectionPoints();
                if plot_iterations
                    figure(1);
                    clf;
                    hold on;
                    plot(cp_initial(1,:), cp_initial(2,:), 'Color',[1 0.5 0.5],'LineWidth',3);
                    plot(spline.cp(1,:), spline.cp(2,:), 'Color',[1 0 0]);

                    gamma = spline.evaluate(t_dense);
                    plot(gamma_initial(1,:), gamma_initial(2,:), 'Color',[0.7 0.7 0.7], 'LineWidth',3);
                    plot(gamma(1,:), gamma(2,:), 'Color',[0 0 0]);

                    scatter(p_infl(1,:), p_infl(2,:), 80, [0 0 1], '.');
                    axis tight equal;
                end
                
                num_infl = length(t_infl);
                if num_infl <= num_infl_initial-2
                    break;
                end
                
                dt_dq = zeros(2*spline.ncp, 1);
                for j=[0 1]
                    cpi = spline.firstControlPointIndexAt(t_infl(ind_infl+j));
                    dofi = (2*cpi-1):(2*cpi+2*spline.degree);
                    dt_infl_dq = spline.inflectionPointSD(t_infl(ind_infl+j));
                    dt_dq(dofi) = dt_dq(dofi) + (1-2*j) * dt_infl_dq(:);
                end
                
                target_change = 0.1;
                options = optimoptions('quadprog','Display','none');
                step = quadprog(BtB, zeros(2*spline.ncp,1),[],[],[dt_dq'; lin_constr_mat],[target_change; zeros(size(lin_constr_rhs))],[],[],[],options);
            
                spline.cp = spline.cp + reshape(step,2,[]);
            end
        end
        
        function [mat, rhs] = buildEndpointConstraints(spline, endpoints_fixed, tangents_fixed)
            if nargin < 3
                tangents_fixed = [false false];
            end
            t_ends = [0, spline.t_max];
            dgamma = spline.evaluateD(t_ends);
            normals = [-dgamma(2,:); dgamma(1,:)];
            normals(:,1) = normals(:,1) / norm(normals(:,1));
            normals(:,2) = normals(:,2) / norm(normals(:,2));
            cpi_ends = spline.firstControlPointIndexAt(t_ends);
            w_ends = spline.evaluateWeights(t_ends);
            dw_ends = spline.evaluateWeightsD(t_ends);
            lin_constr_mat = zeros(6, 2*spline.ncp);
            c1 = cpi_ends(1):cpi_ends(1)+spline.degree;
            lin_constr_mat(1,2*c1-1) = w_ends(:,1);
            lin_constr_mat(2,2*c1) = w_ends(:,1);
            c2 = cpi_ends(2):cpi_ends(2)+spline.degree;
            lin_constr_mat(3,2*c2-1) = w_ends(:,2);
            lin_constr_mat(4,2*c2) = w_ends(:,2);
            
            lin_constr_mat(5,2*c1-1) = dw_ends(:,1)*normals(1,1);
            lin_constr_mat(5,2*c1) = dw_ends(:,1)*normals(2,1);
            lin_constr_mat(6,2*c2-1) = dw_ends(:,2)*normals(1,2);
            lin_constr_mat(6,2*c2) = dw_ends(:,2)*normals(2,2);
            
            p_ends = spline.evaluate(t_ends);
            lin_constr_rhs = [p_ends(:,1); p_ends(:,2); 0; 0];
            
            mat = lin_constr_mat([endpoints_fixed(:); tangents_fixed(:)],:);
            rhs = lin_constr_rhs([endpoints_fixed(:); tangents_fixed(:)]);
        end
        
        function BtB = buildPreconditioningMatrix(spline, num_samples)
            if nargin < 2
                num_samples = 2e2;
            end
            t_dense = linspace(0,spline.t_max,num_samples);
            
            first_cpi = spline.firstControlPointIndexAt(t_dense);
            dofix = 2*first_cpi + (-1:2:2*spline.degree-1)';
            rowix = repmat(1:2:2*length(t_dense)-1, spline.degree+1, 1);
            
            weights = spline.evaluateWeights(t_dense);
            B = zeros(2*length(t_dense), 2*spline.ncp);
            B(sub2ind(size(B), rowix(:), dofix(:))) = weights(:);
            B(sub2ind(size(B), rowix(:)+1, dofix(:)+1)) = weights(:);
            BtB = B' * B;
        end
        
        function ret = approximateDiscreteCurve(points, degree, num_cp, R_target, t_tol_max, t_tol_min, y_fixed, x_fixed, fix_to_zero, use_min_num_points)
            if nargin < 5
                t_tol_max = 0.25;
                t_tol_min = 0.05;
            end
            if nargin < 7
                y_fixed = [true;true];
            end
            if nargin < 8
                x_fixed = [false;false];
            end
            if nargin < 9
                fix_to_zero = true;
            end
            if nargin < 10
                use_min_num_points = false;
            end
            seg_lens = sqrt(sum((points(:,2:end)-points(:,1:end-1)).^2,1));
            spline = SplineCurve(degree, zeros(2, num_cp));
            t_initial = [0 cumsum(seg_lens)];
            t_initial = t_initial * (spline.t_max / t_initial(end));
            cpi = spline.firstControlPointIndexAt(t_initial);
            dofi = cpi + (0:degree)';
            weights = spline.evaluateWeights(t_initial);
            rowi = repmat(1:size(points,2), degree+1,1);
            initial_solve_mat = sparse(rowi,dofi,weights,size(points,2), num_cp);
            
            cpi_ends = spline.firstControlPointIndexAt([0 spline.t_max]);
            weights_ends = spline.evaluateWeights([0 spline.t_max]);
            
            x_constr_mat = zeros(nnz(x_fixed),num_cp);
            x_constr_ind = 1;
            x_constr_rhs = zeros(nnz(x_fixed),1);
            if x_fixed(1)
                x_constr_mat(x_constr_ind,cpi_ends(1) + (0:degree)') = weights_ends(:,1);
                x_constr_rhs(x_constr_ind) = points(1,1);
                x_constr_ind = x_constr_ind + 1;
            end
            if x_fixed(2)
                x_constr_mat(x_constr_ind,cpi_ends(2) + (0:degree)') = weights_ends(:,2);
                x_constr_rhs(x_constr_ind) = points(1,end);
%                 x_constr_ind = y_constr_ind + 1;
            end
            
            y_constr_mat = zeros(nnz(y_fixed),num_cp);
            y_constr_ind = 1;
            y_constr_rhs = zeros(nnz(y_fixed),1);
            if y_fixed(1)
                y_constr_mat(y_constr_ind,cpi_ends(1) + (0:degree)') = weights_ends(:,1);
                y_constr_rhs(y_constr_ind) = points(2,1);
                y_constr_ind = y_constr_ind + 1;
            end
            if y_fixed(2)
                y_constr_mat(y_constr_ind,cpi_ends(2) + (0:degree)') = weights_ends(:,2);
                y_constr_rhs(y_constr_ind) = points(2,end);
                %y_constr_ind = y_constr_ind + 1;
            end
            
            
            options = optimoptions('lsqlin','Display','none');
            if fix_to_zero
                x_constr_rhs(:) = 0;
                y_constr_rhs(:) = 0;
            end
            if use_min_num_points
%                 ndof_x = num_cp - nnz(x_fixed);
                ind_x = round(linspace(1,size(points,2),num_cp));
                initial_solve_mat_x = initial_solve_mat(ind_x,:);
                cpx0 = lsqlin(initial_solve_mat_x, points(1,ind_x)',[],[],x_constr_mat,x_constr_rhs,[],[],[],options)';
                
%                 ndof_y = num_cp - nnz(y_fixed);
                ind_y = round(linspace(1,size(points,2),num_cp));
                initial_solve_mat_y = initial_solve_mat(ind_y,:);
                cpy0 = lsqlin(initial_solve_mat_y, points(2,ind_y)',[],[],y_constr_mat,y_constr_rhs,[],[],[],options)';
            else
                cpx0 = lsqlin(initial_solve_mat, points(1,:)',[],[],x_constr_mat,x_constr_rhs,[],[],[],options)';
                cpy0 = lsqlin(initial_solve_mat, points(2,:)',[],[],y_constr_mat,y_constr_rhs,[],[],[],options)';
            end
            
            cp_initial = [cpx0;cpy0];
            spline.cp = cp_initial;
            
            t_dense = linspace(0,spline.t_max,2e2);
            gamma_initial = spline.evaluate(t_dense);
            
            first_cpi = spline.firstControlPointIndexAt(t_dense);
            dofix = 2*first_cpi + (-1:2:2*spline.degree-1)';
            rowix = repmat(1:2:2*length(t_dense)-1, spline.degree+1, 1);
            
            weights = spline.evaluateWeights(t_dense);
            B = zeros(2*length(t_dense), 2*spline.ncp);
            B(sub2ind(size(B), rowix(:), dofix(:))) = weights(:);
            B(sub2ind(size(B), rowix(:)+1, dofix(:)+1)) = weights(:);
            BtB = B' * B;
            
            lin_constr_mat = zeros(nnz(y_fixed),num_cp*2);
            lin_constr_mat(:,2:2:end) = y_constr_mat;
            
%             figure;
%             hold on;
%             scatter(points(1,:), points(2,:));
% %             plot(cp_initial(1,:), cp_initial(2,:), 'LineWidth', 2, 'Color', [1 0.5 0.5]);
%             plot(gamma_initial(1,:), gamma_initial(2,:), 'LineWidth', 1, 'Color', [1 0 0]);
%             axis tight equal;
            
%             figure;
%             t_dense = linspace(0,spline.t_max,1e3);
%             k = spline.curvature(t_dense);
%             plot(t_dense, k);
%             axis tight;

            

            max_iter = 200;
            p_approx = 100;
            
            tol_infl = [true true];
            [t_infl, p_infl] = spline.findInflectionPoints([t_tol_min, t_tol_min]);
            if any(t_infl < spline.t_max*0.5)
                tol_infl(1) = true;
            end
            if any(t_infl > spline.t_max*0.5)
                tol_infl(2) = true;
            end
            
%             figure(1);
%             clf;
%             hold on;
%             
%             figure(2);
%             clf;
%             hold off;
            
            K_start = [];
            for i=1:max_iter
                gamma = spline.evaluate(t_dense);
                
%                 figure(1);
%                 clf;
%                 hold on;
%                 plot(cp_initial(1,:), cp_initial(2,:), 'Color',[1 0.5 0.5],'LineWidth',3);
%                 plot(spline.cp(1,:), spline.cp(2,:), 'Color',[1 0 0]);
%                 
%                 plot(gamma_initial(1,:), gamma_initial(2,:), 'Color',[0.7 0.7 0.7], 'LineWidth',3);
%                 plot(gamma(1,:), gamma(2,:), 'Color',[0 0 0]);
%                 axis tight equal;
                
                t_tol = [t_tol_min t_tol_min];
                if tol_infl(1)
                    t_tol(1) = t_tol_max;
                end
                if tol_infl(2)
                    t_tol(2) = t_tol_max;
                end
                kappa = spline.curvature(t_dense);
                [t_infl, p_infl] = spline.findInflectionPoints(t_tol);
                
                if length(t_infl) >= 2
%                     plot(p_infl(1,:), p_infl(2,:), 'Color',[0 0 1], 'LineWidth',2);
                    
                    infl_to = p_infl(:,2) - p_infl(:,1);
                    a = [-infl_to(2); infl_to(1)];
                    b = -dot(p_infl(:,1), a);
                    K = (sum(a.*gamma,1)+b) ./ kappa;
                    if K(1) < 0
                        K = -K;
                        a = -a;
                        b = -b;
                    end
                    if any(K <= 0)
                        fprintf('Inflection line crosses curve.');
                    end
                else
                    lpopt = LPStiffnessOptimizer(gamma, kappa);
                    if length(t_infl) == 1
                        scatter(p_infl(1), p_infl(2), 80, [0 0 1], '.');
                        lpopt.gamma_infl = p_infl;
                        [K, b, a] = lpopt.optimizeWithInflections();
                    else
                        [K, b, a] = lpopt.optimizeSimple();
                    end
                end
                
                [objval, K, objgrad, dK] = SplineCurveTools.computeSoftR(spline, p_approx, a, b, t_dense, t_tol);
                if i==1
                    K_start = K;
                    
%                     gamma = spline.evaluate(t_dense);
%                     to = gamma(:,2:end) - gamma(:,1:end-1);
%                     seg_lens = sqrt(sum(to.^2,1));
%                     s = [0 cumsum(seg_lens)];
%                     strip = GeometryGenerator.generateProfileOutlineLoop(s,K_start/min(K_start)*0.01,0,false);
%                     figure('Color','white');
%                     patch('Vertices',strip','Faces',[1:size(strip,2)-1;2:size(strip,2)]');
%                     axis tight equal;
                end
                
%                 figure(2);
%                 clf;
%                 plot(t_dense, K_start/min(K_start), 'LineWidth',3, 'Color',[0.5 0.5 0.5]);
%                 hold on;
%                 plot(t_dense, K/min(K), 'Color',[0 0 0]);
%                 
%                 fprintf('Obj val: %.2f\n', objval);
                c= [b;a];
                c=c/norm(c);
%                 fprintf('c: (%.4f,%.4f,%.4f)\n',c(1),c(2),c(3));
                if objval < R_target
                    break;
                end
                
                target_change_f = -objval * 0.05;
                options = optimoptions('quadprog','Display','none');
                dq = quadprog(BtB, zeros(2*spline.ncp,1),[],[],...
                    [objgrad'; lin_constr_mat],[target_change_f;zeros(nnz(y_fixed),1)],...
                    [],[],[],options);
%                 dq = objgrad / dot(objgrad, objgrad) * target_change_f;
                spline.cp = spline.cp + reshape(dq, 2, []);
                
                tol_infl = [false false];
                if any(t_infl < spline.t_max*0.5)
                    tol_infl(1) = true;
                end
                if any(t_infl > spline.t_max*0.5)
                    tol_infl(2) = true;
                end
            end
            
            min_K = min(K);
            a = a/min_K;
            b = b/min_K;
            K = K/min_K;
            ret = struct('cp_initial', cp_initial, 'cp_final', spline.cp, 'degree', spline.degree, ...
                'R_soft', objval, 'R_hard', max(K), 'K', K, 'a', b, 'b', a, ...
                'gamma_sampled',spline.evaluate(t_dense),'t_sampled',t_dense);
        end
    end
    
    methods(Static, Access=private)
        function [obj, grad] = controlPointObjective(cp0, x)
            grad = x - cp0;
            obj = 0.5 * sum(grad.^2);
        end
        
        function [c,ceq,gc,gceq] = infl_constr(spline,num_infl,x,line_lin,line_const)
            c = [];
            gc = [];
            spline.cp = reshape(x,2,[]);
            t_infl = spline.findInflectionPoints();
            if num_infl ~= length(t_infl)
                ceq = inf(num_infl,1);
                gceq = inf(length(x), num_infl);
                return;
            end
            
            p_infl = spline.evaluate(t_infl);
            ceq = sum(p_infl.*line_lin, 1)' + line_const;
            cpi = spline.firstControlPointIndexAt(t_infl);
            gceq = zeros(length(x), num_infl);
            for i=1:num_infl
                [ds_dq, dp_dq] = spline.inflectionPointSD(t_infl(i));
                grad = sum(dp_dq.*line_lin, 1)';
                gceq((cpi(i)*2-1):(cpi(i)*2-2+length(grad)), i) = grad;
            end
        end
        
        function [objval, K, objgrad, dK] = computeSoftR(spline, p_approx, a, b, t, t_tolerance)
            if nargin < 6
                t_tolerance = [0 0];
            end
            
            objval = [];
            objgrad = [];
                
            gamma = spline.evaluate(t);
            kappa = spline.curvature(t);
            [t_infl, p_infl] = spline.findInflectionPoints(t_tolerance);
            first_cpi = spline.firstControlPointIndexAt(t);
            weights = spline.evaluateWeights(t);
            
            K = (sum(a .* gamma, 1) + b) ./ kappa;
            
            to = gamma(:,2:end)-gamma(:,1:end-1);
            seg_lens = sqrt(sum(to.^2,1));
            w = 0.5 * [seg_lens(1), seg_lens(1:end-1)+seg_lens(2:end), seg_lens(end)];
            
            
            if nargout < 3
                objval = SplineCurveOptimizer.softMaxOverMin(K,w,p_approx);
                return;
            end
            
            dK1 = zeros(2*(spline.degree+1), length(t));
            dK1(1:2:end) = a(1)*weights ./ kappa;
            dK1(2:2:end) = a(2)*weights ./ kappa;
            dkappa = spline.curvatureSD(t);
            dK2 = -K./kappa.*dkappa;
            dK12 = dK1 + dK2; % dK taking into account changes in gamma and kappa (but not a and b)
            
            dK = zeros(2*spline.ncp, length(t));
            rowi = 2*first_cpi-1 + (0:2*spline.degree+1)';
            coli = repmat(1:length(t), 2*(spline.degree+1), 1);
            dK(sub2ind(size(dK), rowi(:), coli(:))) = dK12(:);
            
            
            if length(t_infl)==2
                infl_to = p_infl(:,2) - p_infl(:,1);
                a_cf = [-infl_to(2); infl_to(1)];
                
                [~, dp1_dq] = spline.inflectionPointSD(t_infl(1));
                [~, dp2_dq] = spline.inflectionPointSD(t_infl(2));
                
                infl_first_cpi = spline.firstControlPointIndexAt(t_infl);
                dp1_dq_exp = zeros(2,2*spline.ncp);
                dp1_dq_exp(:,infl_first_cpi(1)*2-1:(infl_first_cpi(1)+spline.degree)*2) = dp1_dq;
                dp2_dq_exp = zeros(2,2*spline.ncp);
                dp2_dq_exp(:,infl_first_cpi(2)*2-1:(infl_first_cpi(2)+spline.degree)*2) = dp2_dq;
                
                dinfl_to = dp2_dq_exp - dp1_dq_exp;
                da = [-dinfl_to(2,:); dinfl_to(1,:)];
                db = -sum(dp1_dq_exp .* a_cf, 1) - sum(p_infl(:,1) .* da, 1);
                if dot(a, a_cf) < 0
                    da = -da;
                    db = -db;
                end
                
                dK3 = (da(1,:)' .* gamma(1,:) + da(2,:)' .*gamma(2,:) + db')./kappa;
                dK = dK + dK3;
            elseif length(t_infl)==1
                [~, dp1_dq] = spline.inflectionPointSD(t_infl(1));
                infl_first_cpi = spline.firstControlPointIndexAt(t_infl);
                dp1_dq_exp = zeros(2,2*spline.ncp);
                dp1_dq_exp(:,infl_first_cpi(1)*2-1:(infl_first_cpi(1)+spline.degree)*2) = dp1_dq;
                db = -sum(dp1_dq_exp .* a, 1);
                dK3 = db'./kappa;
                dK = dK + dK3;
            end
            
            [objval, dobj_dK] = SplineCurveOptimizer.softMaxOverMin(K,w,p_approx);
            objgrad = dK * dobj_dK';
        end
    end
end
