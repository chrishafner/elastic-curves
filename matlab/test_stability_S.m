spline = SplineCurve.import('spline_infl_simple.txt');


% for i=14:25
%     spline.cp(1,i) = 2*spline.cp(1,13) - spline.cp(1,26-i);
%     spline.cp(2,i) = spline.cp(2,26-i);
% end

opt = SplineCurveOptimizer(spline);
% opt.addAllStandardLinearConstraints();
for t = [0 spline.t_max]
    opt.addPositionConstraint(t, [1; 0]);
    opt.addPositionConstraint(t, [0; 1]);
    opt.addTangentConstraint(t, [1; 0]);
end

opt.addInflectionPointLineConstraint();
opt.satisfyConstraintsWithUnderdeterminedNewton(1.e-10);
% opt.addArcLengthConstraint();

% figure(1);
% clf;
% hold on;
% spline.plotCurve([0 0 0]);
% scatter(spline.cp(1,:), spline.cp(2,:));
% axis equal tight;
% 
% figure(2);
% clf;
% hold on;
% [s_mid,detM] = computeStabilityCurve(spline, 20);
% plot(s_mid, detM);

cp_iter = opt.optimizeStability(50, 3., 20, 0.01, 20, 10000,false);
% % opt.optimizeStability(100, 3., 20, 0.03, 20); % for 'spline_infl_simple.txt'
% opt.spline.writeToFile('spline_infl4_opt_rK6.txt');

% spline = SplineCurve.import('spline_infl4_opt_rK6.txt');
% figure(1);
% spline.plotCurve([1 0 0]);
% axis equal tight;
% 
% figure(2);
% [s_mid,detM] = computeStabilityCurve(spline, 20);
% plot(s_mid, detM);


spline.cp = spline.cp * 3 * 1e-2;

t = linspace(0,spline.t_max,1e3);
gamma = spline.evaluate(t);
to = gamma(:,2:end) - gamma(:,1:end-1);
seg_lens = sqrt(sum(to.^2,1));
s = [0 cumsum(seg_lens)];
total_len = sum(seg_lens);

% figure;
% curv = spline.curvature(t);
% plot(s,curv);
% 
% figure;
% plot(gamma(1,:), gamma(2,:));
% axis tight equal;

lpopt = LPStiffnessOptimizer(spline, 1e3);
K = lpopt.optimizeWithInflections();
% figure;
% plot(s,K);
outer_width = 0.1;
strip = GeometryGenerator.generateProfileOutlineLoop(s,0.5 * outer_width * K'/max(K),5.e-3,false);
% figure;
% plot(strip(1,:), strip(2,:));
% axis tight equal;

save('paper-stability-S/S-spline.mat','spline');

function [s_mid,detM] = computeStabilityCurve(spline, tess)
    [elastica, alpha] = spline.convertToElastica(tess);
    lambda = elastica.findOptimalLambda(alpha);
    s_mid = elastica.sMid();
    detM = elastica.computeConstrainedStabilityMatrix(alpha, lambda);
    
    [alpha, lambda] = elastica.optimizeWithNewton(alpha, lambda);
    ret = elastica.computeAll(alpha, lambda);
    % Stability with eigenvalues (bäh!)
    tic;
    H = ret.H;
    Haa = H(3:end, 3:end);
    B = full(H([1 2], 3:end));

    NB = null(B);
    lhs = NB' * Haa * NB;
    rhs = NB' * NB;

    [V,D] = eig(lhs,rhs);
    D = diag(D);
    ineg = find(D<=0);
    neg_vals = D(D<=0);
    fprintf("Neg eigenvalues: %.4f\n", length(neg_vals));
    toc;
    
%     figure;
%     hold on;
%     elastica.draw(alpha);
%     
%     for ind=ineg(:)'
%         eigvec = NB * V(:,ind);
%         eigvec = eigvec/norm(eigvec);
%         alpha_pert = alpha + eigvec * 1;
%         elastica.draw(alpha_pert,[1 0 0]);
%     end
%     axis tight equal;
    
end