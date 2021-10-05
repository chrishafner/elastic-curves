spline = SplineCurve.import('spline_infl_simple.txt');

opt = SplineCurveOptimizer(spline);
% opt.addAllStandardLinearConstraints();
for t = [0 spline.t_max]
    opt.addPositionConstraint(t, [1; 0]);
    opt.addPositionConstraint(t, [0; 1]);
    opt.addTangentConstraint(t, [1; 0]);
end

opt.addInflectionPointLineConstraint();
opt.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);

num_samples = 4e2;
lpopt = LPStiffnessOptimizer(spline, num_samples);
K = lpopt.optimizeWithInflections();

dgamma_ends = spline.evaluateD([0, spline.t_max]);
alpha_ends = atan2(dgamma_ends(2,:), dgamma_ends(1,:));
gamma = spline.evaluate(linspace(0,spline.t_max,4e2));
[alpha, a0, a1] = AbsoluteAngleElastica.computeAngles(gamma, alpha_ends(1), alpha_ends(2));
gamma_to = gamma(:,2:end)-gamma(:,1:end-1);
seg_lens = sqrt(sum(gamma_to.^2,1));
s = [0 cumsum(seg_lens)];
s_mid = 0.5 * (s(1:end-1)+s(2:end));

elastica = AbsoluteAngleElastica(seg_lens', K, a0, a1, gamma(:,1));
    elastica.addEndPointConstraint(gamma(:,end));
[alpha_init, lambda] = elastica.optimizeWithNewton(alpha);

detM = elastica.computeConstrainedStabilityMatrix(alpha_init, lambda);


fig1 = figure();
clf;
elastica.draw(alpha,'k','-',1);
hold on;
axis tight equal;
title('Fig. 9 (left)');

fig2 = figure();
clf;
ax2 = axes();
plot(s_mid, detM, 'k');
hold on;
title('Fig. 9 (right)');

i_constr = 186;
for dir=[-1.1 1]
    elastica.addLinearPointConstraint(i_constr, [0;1], gamma(2,i_constr));

    
    figure(fig1);
    num_steps = 5;
    max_defl = 4.5 * dir;
    alpha = alpha_init;
    lambda = [0;0;0];
    for i=1:num_steps
        elastica.pc_value(end) = max_defl/(num_steps+1)*i;
        [alpha_new, lambda_new] = elastica.optimizeWithNewton(alpha, lambda);
        if mod(i,2)==0
            elastica.draw(alpha_new, 'r', '-', 1);
            axis tight equal;
        end
        alpha = alpha_new;
        lambda = lambda_new;
    end

    elastica.pc_index = elastica.pc_index(1:end-1);
    elastica.pc_coeff = elastica.pc_coeff(:,1:end-1);
    elastica.pc_value = elastica.pc_value(1:end-1);
    [alpha_final, lambda_final] = elastica.optimizeWithNewton(alpha);
    elastica.draw(alpha_final, 'b', '-', 1);
    axis tight equal;
    
    detM = elastica.computeConstrainedStabilityMatrix(alpha_final, lambda_final);
    figure(fig2);
    plot(s_mid, detM, 'r');
    
end

ax2.YLim = [-750 2000];
