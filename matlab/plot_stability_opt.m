spline1 = SplineCurve.import('spline_infl_simple.txt');

opt = SplineCurveOptimizer(spline1);
opt.addAllStandardLinearConstraints();
opt.addInflectionPointLineConstraint();
opt.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);

cp_iter = opt.optimizeStability(50, 2.8, 40, 0.01, 20, 10000, false);

num_iter = size(cp_iter,3);
i_plot = round(linspace(1,num_iter,5));


fig1=figure();
clf;
ax1=axes();
hold on;
title('Fig. 10 (top left)');

fig2=figure();
clf;
ax2=axes();
hold on;
title('Fig. 10 (top right)');

fig3=figure();
clf;
ax3=axes();
hold on;
title('Fig. 10 (bottom left)');

fig4=figure();
clf;
ax4=axes();
hold on;
title('Fig. 10 (bottom right)');

num_samples = 4e2;
for j=1:length(i_plot)
    i = i_plot(j);
    
    spline = SplineCurve(spline1.degree, cp_iter(:,:,i));
    lpopt = LPStiffnessOptimizer(spline, num_samples);
    K = lpopt.optimizeWithInflections();

    dgamma_ends = spline.evaluateD([0, spline.t_max]);
    alpha_ends = atan2(dgamma_ends(2,:), dgamma_ends(1,:));
    gamma = spline.evaluate(linspace(0,spline.t_max,num_samples));
    curv = spline.curvature(linspace(0,spline.t_max,num_samples));
    [alpha, a0, a1] = AbsoluteAngleElastica.computeAngles(gamma, alpha_ends(1), alpha_ends(2));
    gamma_to = gamma(:,2:end)-gamma(:,1:end-1);
    seg_lens = sqrt(sum(gamma_to.^2,1));
    s = [0 cumsum(seg_lens)];
    s_mid = 0.5 * (s(1:end-1)+s(2:end));

    elastica = AbsoluteAngleElastica(seg_lens', K, a0, a1, gamma(:,1));
        elastica.addEndPointConstraint(gamma(:,end));
    [alpha_init, lambda] = elastica.optimizeWithNewton(alpha);

    detM = elastica.computeConstrainedStabilityMatrix(alpha_init, lambda);
    
    figure(fig1);
    plot(gamma(1,:), gamma(2,:), 'k');
    axis tight equal;
    
    figure(fig2);
    plot(s_mid, detM, 'k');
    
    figure(fig3);
    plot(s, curv, 'Color', [1 (length(i_plot)-j)/(length(i_plot))*ones(1,2)]);
    
    figure(fig4);
    plot(s, K, 'Color', [1 (length(i_plot)-j)/(length(i_plot))*ones(1,2)]);

end

ax2.XLim(2) = 13;
ax2.YLim = [-300 200];
ax3.XLim(2) = 13;
ax4.XLim(2) = 13;
ax4.YLim(1) = 0;