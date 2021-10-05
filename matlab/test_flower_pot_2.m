load('flower-pot/curves_R2.5_ncp10_v2.mat', 'curves_opt');
num_curves = length(curves_opt);

special_param = 9.5878;
g = 9.81;
e_gravity = [0; 12*g*special_param];

% Uniform scaling
scale = 6e-2;
num_samples = 500;
gamma_all = zeros(2,num_samples,num_curves);
K_all = zeros(num_curves, num_samples);

% Compute stiffnesses under gravity
% figure;
for ci=1:num_curves
    curves_opt(ci).origin = scale * curves_opt(ci).origin;
    curves_opt(ci).cp_initial = scale * curves_opt(ci).cp_initial;
    curves_opt(ci).cp_final = scale * curves_opt(ci).cp_final;
    curves_opt(ci).gamma_sampled = scale * curves_opt(ci).gamma_sampled;
    curves_opt(ci).a = scale * curves_opt(ci).a;
    curves_opt(ci).b = scale * curves_opt(ci).b;
    
    curve = curves_opt(ci);
    
    spline = SplineCurve(curves_opt(ci).degree, curves_opt(ci).cp_final);
    [t_infl, p_infl] = spline.findInflectionPoints();
    t = linspace(0, spline.t_max, num_samples);
    gamma = spline.evaluate(t);
    curv = spline.curvature(t);
    gamma_to = gamma(:,2:end) - gamma(:,1:end-1);
    seg_lengths = sqrt(sum(gamma_to.^2, 1));
    total_length = sum(seg_lengths);
    s = [0, cumsum(seg_lengths)]';
    eff_lengths = 0.5 * [seg_lengths(1), seg_lengths(1:end-1) + seg_lengths(2:end), seg_lengths(end)];
    
    opt = LPStiffnessOptimizer(gamma, curv);
    opt.gamma_infl = p_infl;
    opt.v_weights = eff_lengths;
    opt.e_gravity = e_gravity;
    
    K_nograv = opt.optimizeWithInflections();
    [K_grav, a_grav, b_grav] = opt.optimizeWithGravity();
    M_opt = max(K_grav);
    [K_fine, a_fine, b_fine] = opt.fineTuneWithGravity(M_opt, 0.01);
    
    gamma_all(:,:,ci) = gamma;
    K_all(ci,:) = K_fine;
    
    fprintf('Strip %u: R=%f\n', ci, max(K_fine)/min(K_fine));
    
%     plot(s,curv);
%     hold on;
    
%     dg = spline.evaluateD([0 spline.t_max]);
%     a0 = atan2(dg(2,1), dg(1,1));
%     a1 = atan2(dg(2,2), dg(1,2));
%     [alpha_init, alpha_start, alpha_end] = AbsoluteAngleElastica.computeAngles(gamma, a0, a1);
%     elastica = AbsoluteAngleElastica(seg_lengths', K_nograv, alpha_start, alpha_end, gamma(:,1));
%     elastica.addEndPointConstraint(gamma(:,end));
%     alpha = elastica.optimizeWithNewton(alpha_init, [0;0]);
%     elastica.e_gravity = e_gravity;
%     alpha_grav = elastica.optimizeWithNewton(alpha_init, [0;0]);
%     
%     elastica.stiffnesses = K_fine;
%     elastica.recomputeDerived();
%     alpha_grav_corrected = elastica.optimizeWithNewton(alpha_init, [0;0]);
%     
%     figure(5);
%     elastica.draw(alpha_init,[0.5 0.5 0.5]);
%     hold on;
%     elastica.draw(alpha,[1 0.5 0.5]);
%     elastica.draw(alpha_grav,[1 0 0]);
%     elastica.draw(alpha_grav_corrected,[0 0 1]);
%     axis equal tight;
end


% Plot all optimized curves
% figure('Color','white');
% for ci=1:num_curves
%     curve = curves_opt(ci);
%     p = curve.origin + curve.plane * gamma_all(:,:,ci);
%     patch('Vertices',p','Faces',[1:size(p,2)-1;2:size(p,2)]');
% end
% axis equal tight vis3d;

% Compute projection of curves onto ground plane and angle bisectors
curve_lines = zeros(3, num_curves);
bisector_lines = zeros(3, num_curves);
for ci=1:num_curves
    o1 = curves_opt(ci).origin;
    d1 = curves_opt(ci).plane(:,1);
    n1 = [-d1(2);d1(1)];
    a1 = -dot(n1,o1([1 2]));
    curve_lines(:,ci) = [n1;a1];
end
for ci=1:num_curves
    ci_next = mod(ci,num_curves)+1;
    na1 = curve_lines(:,ci);
    na2 = curve_lines(:,ci_next);
    na_bis = na1+na2;
    n_bis_norm = norm(na_bis([1 2]));
    bisector_lines(:,ci) = na_bis/n_bis_norm;
end

% Compute max half-width of each curve point (min of normal distances to bisector
% lines)
curve_hw = zeros(num_curves, num_samples);
hw_reduce = 2.e-3;
for ci=1:num_curves
    curve = curves_opt(ci);
    p = curve.origin([1 2]) + curve.plane([1 2],:) * gamma_all(:,:,ci);
    n = curve_lines([1 2],ci);
    bis1 = bisector_lines(:,ci);
    bis2 = bisector_lines(:,mod(ci-2, num_curves)+1);
    d1 = (-bis1(3)-sum(p.*bis1([1 2]),1)) ./ sum(n.*bis1([1 2]), 1);
    d1 = abs(d1) - hw_reduce ./ abs(sum(n.*bis1([1 2]), 1));
    d2 = (-bis2(3)-sum(p.*bis2([1 2]),1)) ./ sum(n.*bis2([1 2]), 1);
    d2 = abs(d2) - hw_reduce ./ abs(sum(n.*bis2([1 2]), 1));
    curve_hw(ci,:) = min([d1; d2],[],1);
end

strip_min_dist = 5e-3;
hw_max = max(curve_hw,[],2)+strip_min_dist*0.5;
strip_offsets = [0; cumsum(hw_max(1:end-1)+hw_max(2:end))];

% Compute, plot, and export glued double-strip design
paper_thickness_1 = 0.2e-3;
paper_thickness_2 = 0.1e-3;
curve_hw_under = zeros(size(curve_hw));

figure;
offset_x = 0.3;
svg_cuts = cell(6*num_curves,1);
slit_length = 2.e-3;
for ci=1:num_curves
    curve_hw_under(ci,:) = GeometryGenerator.computeSubstripHalfWidth(paper_thickness_1,paper_thickness_2,curve_hw(ci,:),K_all(ci,:));
    
    to = gamma_all(:,2:end,ci) - gamma_all(:,1:end-1,ci);
    s = [0, cumsum(sqrt(sum(to.^2,1)))];
    loop_outer = GeometryGenerator.generateProfileOutlineLoop(s,curve_hw(ci,:),6.e-3) + [0; -strip_offsets(ci)];
    loop_inner = GeometryGenerator.generateProfileOutlineLoop(s,curve_hw_under(ci,:),6.e-3) + [offset_x; -strip_offsets(ci)];
    
    patch('Vertices', loop_outer', 'Faces', [1:size(loop_outer,2); 2:size(loop_outer,2) 1]');
    patch('Vertices', loop_inner', 'Faces', [1:size(loop_inner,2); 2:size(loop_inner,2) 1]', 'EdgeColor', 'r');
    
    svg_cuts{6*(ci-1)+1} = loop_outer;
    svg_cuts{6*(ci-1)+2} = loop_inner;
    
    x_slits = [min(loop_outer(1,:)), max(loop_outer(1,:)), min(loop_inner(1,:)), max(loop_inner(1,:))];
    x_dir = [1 -1 1 -1];
    
    for i=1:4
        svg_cuts{6*(ci-1)+2+i} = [x_slits(i) x_slits(i)+x_dir(i)*slit_length; -strip_offsets(ci)*ones(1,2)];
    end
end
axis tight equal;
title('Fig. 19 (right)');

SvgTools.exportCurves('flower-pot/flower-pot-double-strip.svg', svg_cuts, 1e3/0.352778);

% 
% 
% % Compute and plot single-strip simple design
% svg_cuts = cell(num_curves, 1);
% figure;
% for ci=1:num_curves
%     c = (curve_hw(ci,:) ./ K_all(ci,:));
%     hw = K_all(ci,:) * min(c);
%     to = gamma_all(:,2:end,ci) - gamma_all(:,1:end-1,ci);
%     s = [0, cumsum(sqrt(sum(to.^2,1)))];
%     
%     loop = GeometryGenerator.generateProfileOutlineLoop(s,hw,6.e-3) + [0; -(ci-1)*30.e-3];
%     svg_cuts{ci} = loop;
%     
%     patch('Vertices', loop', 'Faces', [1:size(loop,2); 2:size(loop,2) 1]');
% end
% axis tight equal;
% SvgTools.exportCurves('flower-pot/flower-pot-single-strip.svg', svg_cuts, 1e3/0.352778);
% 
% figure('Color','white');
% % Visualize strip widths
% for ci=1:num_curves
%     curve = curves_opt(ci);
%     p = curve.origin + curve.plane * gamma_all(:,:,ci);
%     n = [curve_lines([1 2],ci);0];
%     offset = n.*curve_hw(ci,:);
% 
%     q = [p + offset, p - offset];
%     f = zeros(3, (size(p,2)-1)*2);
%     o = size(p,2);
%     f(1,1:2:end) = 1:o-1;
%     f(2,1:2:end) = 2:o;
%     f(3,1:2:end) = o+2:2*o;
%     f(1,2:2:end) = 1:o-1;
%     f(2,2:2:end) = o+2:2*o;
%     f(3,2:2:end) = o+1:2*o-1;
%     patch('Vertices', q', 'Faces', f', 'EdgeColor','none', 'FaceColor',[0.5 0.5 1]);
%     
%     q2 = [p + offset, fliplr(p - offset)];
%     patch('Vertices', q2', 'Faces', [1:size(q,2); [2:size(q,2) 1]]', 'EdgeColor',[0 0 0]);
% end
% axis tight equal vis3d;
% camlight;
% 
% 
% figure;
% hold on;
% cut_index = 1;
% % Generate perforated strips
% 
% fid = fopen('flower-pot/ring_base.txt', 'w');
% fs = FsTools(fid);
% 
% strips = cell(num_curves,1);
% for ci=1:num_curves
%     curve = curves_opt(ci);
%     K = K_all(ci,:);
%     to = gamma_all(:,2:end,ci) - gamma_all(:,1:end-1,ci);
%     s = [0, cumsum(sqrt(sum(to.^2,1)))];
%     outer_width = 2*curve_hw(ci,:);
%     c_max = min(outer_width ./ K);
%     width_profile = K * c_max;
%     loops = GeometryGenerator.generateMicroPerforationProfile(s,width_profile,outer_width,2,3,0.2,5e-3);
%     if (ci==1)
%         all_cuts = cell(length(loops)*num_curves,1);
%     end
%     loops_offset = cell(size(loops));
%     for i=1:length(loops)
%         loops_offset{i} = (loops{i} + [0;strip_offsets(ci)])*1e2;
% %         if ci < 5
%             plot(loops_offset{i}(1,:), loops_offset{i}(2,:),'k');
% %         end
%         all_cuts{cut_index} = loops_offset{i};
%         cut_index = cut_index + 1;
%     end
%     strips{ci} = loops_offset;
%     
%     % Export single strip
% %     filename_single = sprintf('flower-pot/single/ring%u.svg', ci);
% %     SvgTools.exportCurves(filename_single, loops, 10/0.352778);
%     
%     spline = SplineCurve(curve.degree, curve.cp_final);
%     t_ends = [0 spline.t_max];
%     p_ends = curve.origin + curve.plane * spline.evaluate(t_ends);
%     tangent_ends = curve.plane * spline.evaluateD(t_ends);
%     tangent_ends = tangent_ends ./ sqrt(sum(tangent_ends.^2,1));
%     tangent_ends(:,2) = -tangent_ends(:,2);
%     plane_normal = [curve_lines([1 2], ci); 0];
%     hw_ends = 0.5 * outer_width([1 end]);
%     R = zeros(3,3,2);
%     slot_depth = 6e-3;
%     paper_width = 0.7e-3;
%     for i=1:2
%         R(:,1,i) = plane_normal*(3-2*i);
%         R(:,3,i) = tangent_ends(:,i);
%         R(:,2,i) = cross(R(:,3,i), R(:,1,i));
%         
%         fs.writeCuboid(1e2*[-hw_ends(i) -paper_width/2 -slot_depth], 1e2*[hw_ends(i) paper_width/2 0], 1e2*p_ends(:,i), R(:,:,i));
%     end
% end
% fclose(fid);
% axis equal tight;
% 
% % SvgTools.exportCurves('flower-pot/glass_ring_cuts.svg', all_cuts, 10/0.352778);
% % SvgTools.exportForCricut('flower-pot/flower_pot_cricut.svg', strips, 10/0.352778);
