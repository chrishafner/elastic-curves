% SI units
special_param_paper = 9.5878;
g = 9.81;
num_samples = 3e2;
e_gravity_paper = [0; 12*g*special_param_paper];

curves = struct('gamma',[],'dgamma',[],'kappa',[],'K',[],'vertex_weights',[],'e_gravity',[]);

% Arc from pavilion
ellipse_radii = [75.e-3, 112.5e-3];
num_turns = 1.5;
num_strips = 31;
center_height = 20.e-3;
si = 4;

ellipse_angle = 2*pi*num_turns*(si-1)/(num_strips-1);
ellipse_basis = [cos(ellipse_angle), -sin(ellipse_angle); ...
    sin(ellipse_angle), cos(ellipse_angle)] .* ellipse_radii;
gamma_func = @(t) ellipse_basis * [cos(t);sin(t)] + [0;center_height];
dgamma_func = @(t) ellipse_basis * [-sin(t);cos(t)];
ddgamma_func = @(t) -gamma_func(t);
curv_func = @(t) curv_func_aux(dgamma_func, ddgamma_func, t);

t_samp = linspace(0,2*pi,40);
gamma = gamma_func(t_samp);
i1 = find(gamma(2,1:end)<=0 & gamma(2,[2:end 1])>0);
i2 = find(gamma(2,1:end)>0 & gamma(2,[2:end 1])<=0);
ts = [t_samp(i1), t_samp(i2)];
if ts(2) < ts(1)
    ts(2) = ts(2) + 2*pi;
end
for curve_index=[1 2]
    [ts(curve_index), converged] = NumOpt.newtonsMethod(...
        @(t) sum(gamma_func(t).*[0;1],1), ...
        @(t) sum(dgamma_func(t).*[0;1],1), ...
        ts(curve_index), 1.e-12, 10);
end

t_samp = linspace(ts(1),ts(2),num_samples);
gamma = gamma_func(t_samp);
dgamma = dgamma_func(t_samp);
kappa = curv_func(t_samp);
dgamma_ends = dgamma_func(ts);
a0 = atan2(dgamma_ends(2,1), dgamma_ends(1,1));
a1 = atan2(dgamma_ends(2,2), dgamma_ends(1,2));
[a0, a1] = fixEndpointAngles(gamma, a0, a1);
lpopt = LPStiffnessOptimizer(gamma,kappa);
lpopt.e_gravity = e_gravity_paper;
lpopt.v_weights = computeVertexWeights(gamma);
K = lpopt.optimizeWithGravity();

curve_index=1;
curves(curve_index).gamma = gamma;
curves(curve_index).dgamma = dgamma;
curves(curve_index).kappa = kappa;
curves(curve_index).K = K;
curves(curve_index).e_gravity = e_gravity_paper;
curves(curve_index).vertex_weights = lpopt.v_weights;



% Paper spiral
t0 = 4.5*pi;
t1 = 9.5*pi;

scale = 1.e-3;

l = (1/3)*(t1-t0);
l2 = l*l;
l3 = l2*l;
l4 = l3*l;
tangent = -1;
abc = [l4 l3 l2; 4*l3 3*l2 2*l; 12*l2 6*l 2]\[-tangent*l;-tangent;0];
f = @(x) (x<l) .* (abc(1)*x.*x.*x.*x + abc(2)*x.*x.*x + abc(3)*x.*x + tangent*x);
df = @(x) (x<l) .* (4*abc(1)*x.*x.*x + 3*abc(2)*x.*x + 2*abc(3)*x + tangent);
ddf = @(x) (x<l) .* (12*abc(1)*x.*x + 6*abc(2)*x + 2*abc(3));
x = linspace(0,2*l,1e2);
y = f(x);

g = @(t) scale * [cos(t);sin(t)] .* (t + f(t-t0) - f(t1-t));
dg = @(t) scale * ([-sin(t);cos(t)] .* (t + f(t-t0) - f(t1-t)) + [cos(t);sin(t)] .* (1 + df(t-t0) + df(t1-t)));
ddg = @(t) scale * ([cos(t);sin(t)] .* (-(t + f(t-t0) - f(t1-t)) + (ddf(t-t0) - ddf(t1-t))) + 2 * [-sin(t);cos(t)] .* (1 + df(t-t0) + df(t1-t)));

t = linspace(t0,t1,4e2);
gamma = g(t);
dgamma = dg(t);
ddgamma = ddg(t);
curv = (dgamma(1,:) .* ddgamma(2,:) - dgamma(2,:) .* ddgamma(1,:)) ./ (sum(dgamma.^2, 1) .^ 1.5);


to = gamma(:,2:end) - gamma(:,1:end-1);
seg_lens = sqrt(sum(to.^2,1));
v_weights = 0.5 * [seg_lens(1) seg_lens(1:end-1)+seg_lens(2:end) seg_lens(end)];
total_len = sum(seg_lens);
s = [0 cumsum(seg_lens)];

opt = LPStiffnessOptimizer(gamma, curv);
opt.e_gravity = e_gravity_paper;
opt.v_weights = v_weights;
K = opt.optimizeWithGravity();

curve_index = curve_index + 1;
curves(curve_index).gamma = gamma;
curves(curve_index).dgamma = dgamma;
curves(curve_index).kappa = curv;
curves(curve_index).K = K;
curves(curve_index).e_gravity = e_gravity_paper;
curves(curve_index).vertex_weights = v_weights;


% Flower pot curve
ind_flower = [14 1];
for fi=ind_flower
    load('flower-pot/curves_R2.5_ncp10_v2.mat');
    curve = curves_opt(fi);
    scale = 6e-2;
    spline = SplineCurve(curve.degree, curve.cp_final * scale);
    t = linspace(0,spline.t_max, num_samples);

    curve_index = curve_index + 1;
    curves(curve_index).gamma = spline.evaluate(t);
    curves(curve_index).dgamma = spline.evaluateD(t);
    curves(curve_index).kappa = spline.curvature(t);

    opt = LPStiffnessOptimizer(spline, num_samples);
    opt.e_gravity = e_gravity_paper;

    K_grav = opt.optimizeWithGravity();
    M_opt = max(K_grav);
    K = opt.fineTuneWithGravity(M_opt, 0.01);
    curves(curve_index).K = K;
    curves(curve_index).e_gravity = e_gravity_paper;
    curves(curve_index).vertex_weights = opt.v_weights;
end

% Stability S
curve_index = curve_index + 1;

load('paper-stability-S\S-spline.mat','spline');
t = linspace(0,spline.t_max, num_samples);
curves(curve_index).gamma = spline.evaluate(t);
curves(curve_index).dgamma = spline.evaluateD(t);
curves(curve_index).kappa = spline.curvature(t);

opt = LPStiffnessOptimizer(spline, num_samples);
K = opt.optimizeWithInflections();
curves(curve_index).K = K;
curves(curve_index).e_gravity = [];
curves(curve_index).vertex_weights = opt.v_weights;

save('robustness/curves.mat','curves');

function ret = curv_func_aux(dgf, ddgf, t)
    dg = dgf(t);
    ddg = ddgf(t);
    ret = (dg(1,:) .* ddg(2,:) - dg(2,:) .* ddg(1,:)) ./ (sum(dg.^2, 1) .^ 1.5);
end

function [a0, a1] = fixEndpointAngles(gamma, alpha0, alpha1)
    to = gamma(:,2:end)-gamma(:,1:end-1);
    alpha = atan2(to(2,:), to(1,:))';
    for i=2:length(alpha)
        alpha(i) = fixAngle(alpha(i), alpha(i-1));
    end
    a0 = fixAngle(alpha0, alpha(1));
    a1 = fixAngle(alpha1, alpha(end));
end

function w = computeVertexWeights(gamma)
    to = gamma(:,2:end) - gamma(:,1:end-1);
    seg_lens = sqrt(sum(to.^2,1));
    w = 0.5 * [seg_lens(1) seg_lens(1:end-1)+seg_lens(2:end) seg_lens(end)];
end