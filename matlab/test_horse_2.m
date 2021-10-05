load('horse/horse_cross_sections.mat');
nc = length(all_curves);

plot_iterations = false;

curves_opt = struct('cp',[],'degree',[],'R',[],'K',[],'a',[],'b',[],'t',[],'gamma',[],'origin',[],'plane',[]);
curves_opt(nc).cp = [];
% fig = figure();
for ci=1:nc
% for ci=21
    fprintf('Curve %u:\n', ci);
    curve = all_curves(ci);
    p = curve.p;
    i_conv = convhull(p')';
    if i_conv(2)>i_conv(end-1)
        i_conv = i_conv(end:-1:2);
    else
        i_conv = i_conv(1:end-1);
    end
    
    approx = SplineCurveApproximation(p(:,i_conv),6,4);
    approx.reg_weight = 1.e-5;
    approx.max_iter = 250;
    approx.runOptimization(plot_iterations);
    
    spline = approx.spline;
    [t_infl,p_infl] = spline.findInflectionPoints();
    while ~isempty(t_infl)
        SplineCurveTools.removeBoundaryInflectionPoint(spline, true(4,1), t_infl(1) > 0.5*spline.t_max, 0.1, plot_iterations);
        [t_infl,p_infl] = spline.findInflectionPoints();
    end
    
%     SplineCurveTools.optimizeStiffnessRatio(spline, 2, 0.05, 0.01, true(4,1), plot_iterations, 100, 0.05, false);
    SplineCurveTools.optimizeStiffnessRatio(spline, 2, 0.3, 0.2, true(4,1), plot_iterations, 100, 0.05, false);
    
    ns = 2e2;
    lpopt = LPStiffnessOptimizer(spline,ns);
    [K, a, b] = lpopt.optimizeWithInflections();
    
    curves_opt(ci).cp = spline.cp;
    curves_opt(ci).degree = spline.degree;
    curves_opt(ci).R = max(K)/min(K);
    curves_opt(ci).K = K;
    curves_opt(ci).a = a;
    curves_opt(ci).b = b;
    curves_opt(ci).t = linspace(0,spline.t_max,ns);
    curves_opt(ci).gamma = spline.evaluate(curves_opt(ci).t);
    curves_opt(ci).origin = all_curves(ci).origin;
    curves_opt(ci).plane = all_curves(ci).plane;
    
%     if ~isempty(t_infl)
%     figure(fig);
%     clf;
%     hold on;
%     patch('Vertices',p','Faces',[1:size(p,2)-1;2:size(p,2)]','EdgeColor',0.5*ones(3,1));
%     patch('Vertices',p','Faces',[i_conv(1:end-1); i_conv(2:end)]','Marker','.','MarkerSize',10);

    t = linspace(0,spline.t_max,1e3);
    gamma = spline.evaluate(t);
%     patch('Vertices',gamma','Faces',[1:size(gamma,2)-1;2:size(gamma,2)]','EdgeColor','r');
%     patch('Vertices',p_infl','Faces',(1:size(p_infl,2))','Marker','.','MarkerSize',30,'MarkerEdgeColor','r');

%     axis equal tight;
%     end
end

save('horse/horse_curves.mat', 'curves_opt');