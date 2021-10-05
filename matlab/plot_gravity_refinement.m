load('flower-pot/curves_R2.5_ncp10_v2.mat', 'curves_opt');

special_param = 9.5878;
g = 9.81;
e_gravity = [0; 12*g*special_param];
scale = 6e-2;

figure;

i_sim = [1 7 15];
for j=1:length(i_sim)
    ax = subplot(1,3,j);
    
    i = i_sim(j);
    c = curves_opt(i);
    cp = scale * c.cp_final;
    spline = SplineCurve(c.degree, cp);
    
    t = linspace(0,spline.t_max,3e2);
    g = spline.evaluate(t);
    curv = spline.curvature(t);
    gamma_to = g(:,2:end) - g(:,1:end-1);
    seg_lengths = sqrt(sum(gamma_to.^2, 1));
    s = [0 cumsum(seg_lengths)];
    eff_lengths = 0.5 * [seg_lengths(1), seg_lengths(1:end-1) + seg_lengths(2:end), seg_lengths(end)];

    opt = LPStiffnessOptimizer(g, curv);

    opt.v_weights = eff_lengths;
    opt.e_gravity = e_gravity;

    [K_grav, a_grav, b_grav] = opt.optimizeWithGravity();
    M_opt = max(K_grav);
    [K, a_fine, b_fine] = opt.fineTuneWithGravity(M_opt, 0.05);
    
    plot(s, K_grav);
    hold on;
    plot(s, K);
    ax.YLim = [0 3.5];
end
sgtitle('Fig. 8');