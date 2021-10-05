classdef SplineCurveApproximation < handle
    properties
        points
        weights
        t_closest
        p_closest
        ncp
        degree
        
        reg_weight = 1.e-4
        max_iter = 150
        step_size_tol = 1.e-10
        
        spline
        lambda
        constr_mat
        constr_rhs
        constr_proj_mat
        
    end
    
    methods
        function obj = SplineCurveApproximation(points, ncp, degree)
            obj.points = points;
            obj.ncp = ncp;
            obj.degree = degree;
            
            obj.initializeSpline();
%             obj.runOptimization();
        end
        
        function iterations = runOptimization(obj, plot_iterations)
            if nargin < 2
                plot_iterations = false;
            end
            
            bb = max(obj.points,[],2) - min(obj.points,[],2);
            max_cp_step = 0.1*max(bb);
            
            iterations = CellArray();
            
            min_step_size = inf;
            for i=1:obj.max_iter
                if plot_iterations
                figure(1);
                    clf;
                    scatter(obj.points(1,:), obj.points(2,:));
                    hold on;
                    t = linspace(0,obj.spline.t_max,1e3);
                    g = obj.spline.evaluate(t);
                    plot(g(1,:), g(2,:));
                    to = reshape([obj.points; obj.p_closest],2,[]);
                    patch('Vertices',to','Faces',reshape(1:2*size(obj.points,2),2,[])');
                    patch('Vertices',obj.spline.cp','Faces',[1:obj.ncp-1;2:obj.ncp]');
                    axis equal tight;
                end
                
                [f, grad_f, L, grad_L, hess_L] = obj.computeLagrangian();
                
                step = -hess_L\grad_L;
                step = step(1:2*obj.ncp);
                step_size = max(abs(step));
                min_step_size = min(min_step_size, step_size);
                
                iterations.append(obj.spline.cp);
                
                if step_size < obj.step_size_tol
                    break;
                end
                if step_size > max_cp_step
                    step = step / step_size * max_cp_step;
                end
                q_step = reshape(step,2,[]);
                obj.spline.cp = obj.spline.cp + q_step;
                
                obj.computeClosestPoints(obj.t_closest);
            end
            fprintf('Stopped after %u steps\n', i);
            fprintf('Min step size: %e\n', min_step_size);
            
            iterations = iterations.get();
        end
    end
    
    methods(Access=private)
        function initializeSpline(obj)
            seg_lens = sqrt(sum((obj.points(:,2:end)-obj.points(:,1:end-1)).^2,1));
            dummy_spline = SplineCurve(obj.degree, zeros(2, obj.ncp));
            t_initial = [0 cumsum(seg_lens)];
            t_initial = t_initial * (dummy_spline.t_max / t_initial(end));
            discrete_curve = PPHelper.makePiecewiseLinearVector(t_initial,obj.points);
            t_dif = t_initial(2:end) - t_initial(1:end-1);
            obj.weights = 0.5 * [t_dif(1), t_dif(1:end-1)+t_dif(2:end), t_dif(end)];
            
            t_samples = linspace(0, dummy_spline.t_max, obj.ncp*2);
            
            p_fit = ppval(discrete_curve, t_samples);
            cpi = dummy_spline.firstControlPointIndexAt(t_samples);
            dofi = cpi + (0:obj.degree)';
            w = dummy_spline.evaluateWeights(t_samples);
            rowi = repmat(1:length(t_samples), obj.degree+1,1);
            initial_solve_mat = sparse(rowi,dofi,w,length(t_samples), obj.ncp);
            
%             cp0 = (initial_solve_mat \ p_fit')';
            
            obj.constr_rhs = [obj.points(:,1); obj.points(:,end)];
            obj.constr_mat = zeros(4, 2*obj.ncp);
            obj.constr_mat(1, (2*cpi(1)-1):2:(2*obj.degree+2*cpi(1)-1)) = w(:,1);
            obj.constr_mat(2, (2*cpi(1)):2:(2*obj.degree+2*cpi(1))) = w(:,1);
            obj.constr_mat(3, (2*cpi(end)-1):2:(2*obj.degree+2*cpi(end)-1)) = w(:,end);
            obj.constr_mat(4, (2*cpi(end)):2:(2*obj.degree+2*cpi(end))) = w(:,end);
            obj.lambda = zeros(4,1);
            
            obj.constr_proj_mat = eye(2*obj.ncp) - obj.constr_mat' * ((obj.constr_mat * obj.constr_mat') \ obj.constr_mat);
            
            options = optimoptions('lsqlin','Display','none');
            cp0x = lsqlin(initial_solve_mat, p_fit(1,:)',[],[],obj.constr_mat([1 3],1:2:end),obj.constr_rhs([1 3]),[],[],[],options);
            cp0y = lsqlin(initial_solve_mat, p_fit(2,:)',[],[],obj.constr_mat([2 4],2:2:end),obj.constr_rhs([2 4]),[],[],[],options);
            cp0 = [cp0x, cp0y]';
            obj.spline = SplineCurve(obj.degree, cp0);
            obj.computeClosestPoints(t_initial);
            
%             H = full(H);
%             del = 1.e-6;
%             fp = zeros(2*obj.ncp,2);
%             hp = zeros(2*obj.ncp,2*obj.ncp,2);
%             for ind=1:2*obj.ncp
%                 for i=1:2
%                     obj.spline.cp = cp0;
%                     obj.spline.cp(ind) = cp0(ind) + (3-2*i)*del;
%                     obj.computeClosestPoints(obj.t_closest);
%                     [a, gfp] = obj.computeObjective();
%                     fp(ind,i) = a;
%                     hp(:,ind,i) = gfp;
%                 end
%             end
%             grad_f_num = (fp(:,1) - fp(:,2)) / (2*del);
%             hess_f_num = (hp(:,:,1) - hp(:,:,2)) / (2*del);
            
%             p_dense = obj.spline.evaluate(linspace(0,obj.spline.t_max,1e3));
%             figure;
%             scatter(obj.points(1,:), obj.points(2,:));
%             hold on;
%             plot(p_dense(1,:), p_dense(2,:));
%             axis tight equal;
        end
        
        
        function computeClosestPoints(obj, t_guess)
            iter_max = 10;
            f_tol = 1.e-12;
            for i=1:iter_max
                g = obj.spline.evaluate(t_guess);
                dg = obj.spline.evaluateD(t_guess);
                to = obj.points - g;
                f = sum(to.*dg, 1);
                if max(abs(f)) < f_tol || i==iter_max
                    break;
                end
                
                ddg = obj.spline.evaluateD2(t_guess);
                df = sum(to.*ddg, 1) - sum(dg.*dg,1);
                t_guess = t_guess - f ./ df;
            end
            obj.t_closest = t_guess;
            obj.p_closest = g;
        end
        
        function [f, grad_f, hess_ri, hess_ci, hess_val] = computeObjective(obj)
            to = obj.p_closest - obj.points;
            f = 0.5 * sum(sum(to.*to, 1).*obj.weights);
            
            if nargout > 1
                cpi = obj.spline.firstControlPointIndexAt(obj.t_closest);
                w = obj.spline.evaluateWeights(obj.t_closest);
                grad = repmat(obj.weights .* to, obj.degree+1, 1) .* kron(w,[1;1]);
                dofi = (cpi-1)*2 + (1:2*(obj.degree+1))';
                grad_f = accumarray(dofi(:), grad(:));
                
                if nargout > 2
                    dgamma = obj.spline.evaluateD(obj.t_closest);
                    ddgamma = obj.spline.evaluateD2(obj.t_closest);
                    dgamma_sqnorm = sum(dgamma.*dgamma,1);
                    fac = dgamma_sqnorm + sum(to .* ddgamma,1);
                    dw = obj.spline.evaluateWeightsD(obj.t_closest);
                    numer1 = repmat(dgamma, obj.degree+1, 1) .* kron(w,[1;1]);
                    numer2 = repmat(to, obj.degree+1, 1) .* kron(dw,[1;1]);
                    dt_dq = (numer1 + numer2) ./ (-fac);
                    
                    w2 = permute(w, [1 3 2]) .* permute(w, [3 1 2]);
                    term1 = zeros(2*(obj.degree+1), 2*(obj.degree+1), size(obj.points,2));
                    term1(1:2:end,1:2:end,:) = w2;
                    term1(2:2:end,2:2:end,:) = w2;
                    term2 = (permute(dt_dq, [1 3 2]) .* permute(dt_dq, [3 1 2])) .* permute(fac, [3 1 2]);
                    term3 = permute(numer1, [1 3 2]) .* permute(dt_dq, [3 1 2]);
                    term3 = term3 + permute(term3, [2 1 3]);
                    term4 = permute(numer2, [1 3 2]) .* permute(dt_dq, [3 1 2]);
                    term4 = term4 + permute(term4, [2 1 3]);
                    all_terms = permute(obj.weights, [3 1 2]) .* (term1 + term2 + term3 + term4);
                    rowi = repmat(permute(dofi,[1 3 2]), 1, 2*(obj.degree+1), 1);
                    coli = permute(rowi,[2 1 3]);
                    
                    hess_ri = rowi(:);
                    hess_ci = coli(:);
                    hess_val = all_terms(:);
                    
%                     hess_f = sparse(rowi(:), coli(:), all_terms(:), 2*obj.ncp, 2*obj.ncp);
                end
            end
        end
        
        function [f, grad_f, hess_ri, hess_ci, hess_val] = computeRegObjective(obj)
            t = linspace(0,obj.spline.t_max,200);
            t_dif = t(2)-t(1);
            dg = obj.spline.evaluateD(t);
            dg_norm = sqrt(sum(dg.*dg,1));
            w_int = dg_norm .* t_dif;
            
            ddg = obj.spline.evaluateD2(t);
            ddw = obj.spline.evaluateWeightsD2(t);
            cpi = obj.spline.firstControlPointIndexAt(t);
            f = obj.reg_weight * 0.5 * sum(sum(ddg.*ddg, 1) .* w_int);
            
            grad = repmat(ddg .* w_int, obj.degree+1, 1) .* kron(ddw,[1;1]);
            dofi = (cpi-1)*2 + (1:2*(obj.degree+1))';
            grad_f = obj.reg_weight * accumarray(dofi(:), grad(:));
            
            rowi = repmat(permute(dofi,[1 3 2]), 1, 2*(obj.degree+1), 1);
            coli = permute(rowi,[2 1 3]);
            vals = permute(ddw, [1 3 2]) .* permute(ddw.* w_int, [3 1 2]);
            vals_ext = zeros(2*(obj.degree+1), 2*(obj.degree+1), length(t));
            vals_ext(1:2:end,1:2:end,:) = vals;
            vals_ext(2:2:end,2:2:end,:) = vals;
            
            hess_ri = rowi(:);
            hess_ci = coli(:);
            hess_val = obj.reg_weight * vals_ext(:);
        end
        
        function [f, grad_f, L, grad_L, hess_L] = computeLagrangian(obj)
            [f, grad_f, hess_ri, hess_ci, hess_val] = obj.computeObjective();
            [r, grad_r, r_hess_ri, r_hess_ci, r_hess_val] = obj.computeRegObjective();
            f = f + r;
            grad_f = grad_f + grad_r;
            
            q = obj.spline.cp(:);
            g = obj.constr_mat * q - obj.constr_rhs;
            L = f + dot(obj.lambda, g);
            grad_L = [grad_f + obj.constr_mat' * obj.lambda; g];
            
            cri = repmat((2*obj.ncp+1:2*obj.ncp+4)', 2*obj.ncp, 1);
            cci = kron((1:2*obj.ncp)', ones(4,1));
            hess_L = sparse([hess_ri; cri; cci; r_hess_ri], ...
                [hess_ci; cci; cri; r_hess_ci], ...
                [hess_val; obj.constr_mat(:); obj.constr_mat(:); r_hess_val], ...
                2*obj.ncp+4, 2*obj.ncp+4);
        end
        
    end
end