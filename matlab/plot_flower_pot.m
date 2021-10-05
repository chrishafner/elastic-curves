mesh = HalfEdgeMesh.importOBJ('flower-pot/flower-pot.obj');
load('flower-pot/curves_R2.5_ncp10_v2.mat', 'curves_opt');
num_curves = length(curves_opt);

x_offset = 1.e-3;
for i=1:mesh.nv
    mesh.v_traits(i).position = mesh.v_traits(i).position + x_offset;
end
v_pos_all = mesh.allVertices().position();

fhe = mesh.allFaces().halfedge();
fv1 = fhe.from().ind;
fv2 = fhe.to().ind;
fv3 = fhe.next().to().ind;
fv4 = fhe.prev().from().ind;
fv = [fv1;fv2;fv3;fv4];
figure('Color','white');
patch('Vertices', v_pos_all', 'Faces', fv', 'EdgeColor','none','FaceColor',[1 1 1],'FaceLighting','gouraud',...
    'SpecularStrength',0.3,'FaceAlpha',0.3);
axis tight equal off vis3d;
light('Position',[0 0 -1],'Color',[0.3 0.3 0.3])
view(-23,83);
camlight;
view(-128,9);
hold on;

special_param = 9.5878;
g = 9.81;
e_gravity = [0; 12*g*special_param];
scale = 6e-2;

R_vals = zeros(2, num_curves);
for i=1:length(curves_opt)
    c = curves_opt(i);
    cp = cat(3,c.cp_preprocessed,c.cp_final);
    colors = [167 201 167; 52 153 52]/255;
    line_widths = [2 1];
    for j=1:2
        spline = SplineCurve(c.degree, cp(:,:,j));
        [t_infl, p_infl] = spline.findInflectionPoints();
        t = linspace(0,spline.t_max,3e2);
        g = spline.evaluate(t);
        curv = spline.curvature(t);
        g3 = c.origin + c.plane * g;
        plot3(g3(1,:), g3(2,:),g3(3,:),'Color',colors(j,:),'LineWidth',line_widths(j));
        
        g = scale * g;
        curv = curv / scale;
        p_infl = p_infl * scale;
        gamma_to = g(:,2:end) - g(:,1:end-1);
        seg_lengths = sqrt(sum(gamma_to.^2, 1));
        eff_lengths = 0.5 * [seg_lengths(1), seg_lengths(1:end-1) + seg_lengths(2:end), seg_lengths(end)];

        opt = LPStiffnessOptimizer(g, curv);
        
        opt.v_weights = eff_lengths;
        opt.e_gravity = e_gravity;

%         if j==1
%             opt.gamma_infl = p_infl;
%             K = opt.optimizeWithInflections();
%         else
            [K_grav, a_grav, b_grav] = opt.optimizeWithGravity();
            M_opt = max(K_grav);
            [K, a_fine, b_fine] = opt.fineTuneWithGravity(M_opt, 0.01);
%         end
        R_vals(j,i) = max(K)/min(K);
    end
end
title('Fig. 19 (left)');

figure;
ax = axes();
scatter(0:num_curves-1, R_vals(1,:),'Marker','^','MarkerFaceColor','flat','MarkerEdgeColor','none');
hold on;
scatter(0:num_curves-1, R_vals(2,:),'Marker','s','MarkerFaceColor','flat','MarkerEdgeColor','none');
ax.YScale = 'log';
ax.YLim = [1 300];
ax.XLim = [-1 num_curves];
title('Fig. 19 (center)');
