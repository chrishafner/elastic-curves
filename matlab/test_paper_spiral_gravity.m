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
[K,a,b] = opt.optimizeSimple();

special_param = 9.5878;
g = 9.81;
e_gravity = [0; 12*g*special_param];
opt.e_gravity = e_gravity;
opt.v_weights = v_weights;
K_grav = opt.optimizeWithGravity();


a0 = atan2(dgamma(2,1), dgamma(1,1));
a1 = atan2(dgamma(2,end), dgamma(1,end));
[alpha, alpha_start, alpha_end] = AbsoluteAngleElastica.computeAngles(gamma, a0, a1);

elastica = AbsoluteAngleElastica(seg_lens', K, alpha_start, alpha_end, gamma(:,1));
elastica.addEndPointConstraint(gamma(:,end));
elastica.e_gravity = e_gravity;
[alpha_final, lambda_final] = elastica.optimizeWithNewton(alpha, [0;0]);

colors = [60,135,60; 94,46,117]/255;

figure;
plot(gamma(1,:), gamma(2,:));
hold on;
elastica.draw(alpha,colors(1,:),'-',2);
elastica.draw(alpha_final,colors(2,:),'-',2);
axis tight equal;
title('Fig. 7 (top left)');

outer_width = 5.e-2; % 5 cm
Ks = [K_grav K]';

titles = {'Fig. 7 (bottom right)','Fig. 7 (bottom left)'};
for ki=[1 2]
    width_profile = outer_width * (Ks(ki,:) / max(Ks(ki,:)));
    
    p = GeometryGenerator.generateProfileOutlineLoop(s*1e2,width_profile*0.5e2,0.5,false);
    p(1,:) = -p(1,:);
    figure;
    patch('Vertices', p', 'Faces', [1:size(p,2)-1; 2:size(p,2)]', 'EdgeColor',colors(ki,:),'LineWidth',2);
    axis equal tight;
    title(titles{ki});
    SvgTools.exportForCricut(sprintf('paper-spiral/paper-spiral-cricut-varwidth%u.svg', ki), {{p}}, 10/0.352778);
end

