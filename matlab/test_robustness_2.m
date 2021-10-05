load('robustness/curves.mat');
num_curves = length(curves);

curvature_scale = [-1.5e-3, -0.8e-4 -1e-3 -1e-3 -0.7e-3];
num_curvature_lines = [30, 100, 30, 30, 50];

special_param_paper = 9.5878;
g = 9.81;
num_samples = 3e2;
e_gravity_paper = [0; 12*g*special_param_paper];

figure('Color','white');
for ci=1:num_curves
    curve = curves(ci);
    
    % Plot curve and curvature
    subplot(num_curves, 6, (ci-1)*6 + 1);
    plot(curve.gamma(1,:), curve.gamma(2,:), 'Color','k','LineWidth',2);
    hold on;
    
    
    normals = [-curve.dgamma(2,:); curve.dgamma(1,:)];
    normals = normals ./ sqrt(sum(normals.^2,1));
    kappa_curve = curve.gamma + normals .* curve.kappa * curvature_scale(ci);
    plot(kappa_curve(1,:), kappa_curve(2,:), 'Color','k','LineWidth',0.5);
%     patch('Vertices',kappa_curve','Faces',[1:size(curve.gamma,2)-1;2:size(curve.gamma,2)]','EdgeColor','k','LineWidth',0.5);
    
    ncl = num_curvature_lines(ci);
    i_curv = round(linspace(1,size(curve.gamma,2),ncl));
    patch('Vertices',[curve.gamma(:,i_curv), kappa_curve(:,i_curv)]', ...
        'Faces',[1:ncl; ncl+1:2*ncl]','LineWidth',0.5);
    axis tight equal off;
    
    % Plot stiffness
    ax = subplot(num_curves, 6, (ci-1)*6 + 2);
    hold on;
    to = curve.gamma(:,2:end) - curve.gamma(:,1:end-1);
    lens = sqrt(sum(to.^2,1));
    s = [0 cumsum(lens)];
    
    plot(s,curve.K);
    axis([0 s(end) 0 max(curve.K)+1.e-2]);
    ax.Box = 'off';
    ax.YGrid = 'on';
    if max(curve.K)<2
        ax.YTick = 0:0.5:max(curve.K);
    else
        ax.YTick = 0:1:max(curve.K);
    end
    
    % Plot curve with stiffness deviations
    ax = subplot(num_curves, 6, (ci-1)*6 + 3);
    hold on;
    a0 = atan2(curve.dgamma(2,1), curve.dgamma(1,1));
    a1 = atan2(curve.dgamma(2,end), curve.dgamma(1,end));
    alpha_init = [a0, atan2(to(2,:), to(1,:)), a1];
    for i=2:length(alpha_init)
        alpha_init(i) = fixAngle(alpha_init(i), alpha_init(i-1));
    end
    a1 = alpha_init(end);
    alpha_init = alpha_init(2:end-1)';
    
    stiffness_mod = 2.^linspace(-1,1,5);
    stiffness_cols = [1 0 0; 1 0.5 0.5; 0 0 0; 0.5 0.5 1; 0 0 1];
    
    for i=1:5
        K = (curve.K-1)*stiffness_mod(i) + 1;
        elastica = AbsoluteAngleElastica(lens', K, a0, a1, curve.gamma(:,1));
        elastica.addEndPointConstraint(curve.gamma(:,end));
        elastica.e_gravity = curve.e_gravity;
        [alpha, ~, converged] = elastica.optimizeWithNewton(alpha_init, [0;0]);
        if ~converged
           K = (curve.K-1)*(stiffness_mod(i)*0.8) + 1;
            elastica = AbsoluteAngleElastica(lens', K, a0, a1, curve.gamma(:,1));
            elastica.addEndPointConstraint(curve.gamma(:,end));
            elastica.e_gravity = curve.e_gravity;
            [alpha, ~, converged] = elastica.optimizeWithNewton(alpha_init, [0;0]);
        end
        p = elastica.computePoints(alpha);
        plot(p(1,:), p(2,:), 'Color',stiffness_cols(i,:));
    end
    axis tight equal off;
    
    % Plot with wrong gravity
    ax = subplot(num_curves, 6, (ci-1)*6 + 4);
    hold on;
    gravity_mod = 2.^linspace(-1,1,5);
    for i=1:5
        lpopt = LPStiffnessOptimizer(curve.gamma, curve.kappa);
        lpopt.v_weights = curve.vertex_weights;
        if isempty(curve.e_gravity)
            K = curve.K;
        else
            lpopt.e_gravity = curve.e_gravity * gravity_mod(i);
            K = lpopt.optimizeWithGravity();
        end

        elastica = AbsoluteAngleElastica(lens', K, a0, a1, curve.gamma(:,1));
        elastica.addEndPointConstraint(curve.gamma(:,end));
        if isempty(curve.e_gravity)
            elastica.e_gravity = 0.1 * e_gravity_paper * log(gravity_mod(i))/log(2);
        else
            elastica.e_gravity = curve.e_gravity;
        end
        alpha = elastica.optimizeWithNewton(alpha_init, [0;0]);
        p = elastica.computePoints(alpha);
        plot(p(1,:), p(2,:), 'Color',stiffness_cols(i,:));
    end
    axis tight equal off;
    
    % Plot with wrong tangents
    ax = subplot(num_curves, 6, (ci-1)*6 + 5);
    hold on;
    angle_mod = linspace(-10,10,5)/180*pi;
    if curve.gamma(1,1) > curve.gamma(1,end)
        angle_mod = -angle_mod;
    end
    for i=1:5
        elastica = AbsoluteAngleElastica(lens', curve.K, a0+angle_mod(i), a1-angle_mod(i), curve.gamma(:,1));
        elastica.addEndPointConstraint(curve.gamma(:,end));
        elastica.e_gravity = curve.e_gravity;
        [alpha, ~, converged] = elastica.optimizeWithNewton(alpha_init, [0;0],10,1.e-8);
        if ~converged
            elastica = AbsoluteAngleElastica(lens', curve.K, a0+angle_mod(i)*0.75, a1-angle_mod(i)*0.75, curve.gamma(:,1));
            elastica.addEndPointConstraint(curve.gamma(:,end));
            elastica.e_gravity = curve.e_gravity;
            [alpha, ~, converged] = elastica.optimizeWithNewton(alpha_init, [0;0],10,1.e-8);
        end
        p = elastica.computePoints(alpha);
        plot(p(1,:), p(2,:), 'Color',stiffness_cols(i,:));
    end
    axis tight equal off;
    
    % Plot with plasticity
    plast_cols = flipud([0 0 0; 1 0.75 0.75; 1 0.5 0.5; 1 0.25 0.25; 1 0 0]);
    ax = subplot(num_curves, 6, (ci-1)*6 + 6);
    hold on;
    kappa_max = max(abs(curve.kappa));
    curv_mod = linspace(0.5,1,5);
    for i=1:5
        kappa_natural = zeros(size(curve.K));
        kappa_plast = kappa_max * curv_mod(i);
        kappa_natural(curve.kappa > kappa_plast) = curve.kappa(curve.kappa > kappa_plast) - kappa_plast;
        kappa_natural(curve.kappa < -kappa_plast) = curve.kappa(curve.kappa < -kappa_plast) + kappa_plast;
        
        elastica = AbsoluteAngleElastica(lens', curve.K, a0, a1, curve.gamma(:,1));
        elastica.addEndPointConstraint(curve.gamma(:,end));
        elastica.e_gravity = curve.e_gravity;
        
        elastica.kappa_natural = kappa_natural;
        alpha = elastica.optimizeWithNewton(alpha_init, [0;0]);
        p = elastica.computePoints(alpha);
        plot(p(1,:), p(2,:), 'Color',stiffness_cols(i,:));
    end
    axis tight equal off;
end
sgtitle('Fig. 15');