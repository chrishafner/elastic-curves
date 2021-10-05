mesh = HalfEdgeMesh.importOBJ('teaser/teaser_sixth.obj');
render_mesh = HalfEdgeMesh.importOBJ('teaser/teaser-render.obj');
load('teaser/teaser_curves.mat','curves_opt');
load('teaser/teaser_half_width.mat','hw_all','between_planes');

scale = 0.03;
num_strips = length(curves_opt);
num_samples = size(curves_opt(1).gamma,2);

F = GeometryGenerator.generateStripFaceList(num_samples);

fig1 = figure();
fig1.Color='w';
clf;
hold on;
title('Fig. 20 (center)');
fig2=figure();
fig2.Color='w';
clf;
hold on;
title('Fig. 20 (right)');
light_gray = [216 216 216]/255;
gray = [43 48 43]/255;
green = [60 135 60]/255;
purple = [94 46 117]/255;
light_green = [177 207 177]/255;
    
for si=1:num_strips
    c = curves_opt(si);
    origin = scale * c.origin;
    gamma = scale * c.gamma;
    hw = hw_all(si,:);
    p = origin + c.basis(:,[1 2]) * gamma;
    normal = c.basis(:,3);
    
    p_strip = [p - normal * hw, p + normal * hw];
    
    figure(fig1);
    patch('Vertices',p_strip','Faces',F,'EdgeColor','none','FaceColor',green,'SpecularStrength',0.3,'DiffuseStrength',1.);
    i_bdry = [1:num_samples 2*num_samples:-1:num_samples+1];
    patch('Vertices',p_strip','Faces',[i_bdry; i_bdry([2:end 1])]', 'EdgeColor',gray);
    patch('Vertices', p_strip([1 2],:)', 'Faces',F, 'EdgeColor','none','FaceColor',light_gray,...
        'FaceLighting','none');
    
    figure(fig2);
    hw_K = min(hw ./ c.K') .* c.K';
    pK_strip = [p - normal * hw_K, p + normal * hw_K];
    patch('Vertices',pK_strip','Faces',F,'EdgeColor','none','FaceColor',purple,'SpecularStrength',0.3,'DiffuseStrength',1.);
    patch('Vertices',pK_strip','Faces',[i_bdry; i_bdry([2:end 1])]', 'EdgeColor',gray);
    patch('Vertices', pK_strip([1 2],:)', 'Faces',F, 'EdgeColor','none','FaceColor',light_gray,...
        'FaceLighting','none');
    
    
%     gamma_to = gamma(:,2:end)-gamma(:,1:end-1);
%     seg_lens = sqrt(sum(gamma_to.^2,1));
%     
%     spline = SplineCurve(c.degree, c.cp*scale);
%     dgamma_ends = spline.evaluateD([0,spline.t_max]);
%     angles_ends = atan2(dgamma_ends(2,:), dgamma_ends(1,:));
%     
%     [alpha, alpha_start, alpha_end] = AbsoluteAngleElastica.computeAngles(gamma, angles_ends(1), angles_ends(2));
%     
%     
%     elastica = AbsoluteAngleElastica(seg_lens', hw_all(si,:)', alpha_start, alpha_end, gamma(:,1));
%     elastica.addEndPointConstraint(gamma(:,end));
%     alpha_final = elastica.optimizeWithNewton(alpha);
%     p = elastica.computePoints(alpha_final);
%     p3d = origin + c.basis(:,[1 2])*p;
%     patch('Vertices',p3d','Faces',[1:size(p3d,2)-1;2:size(p3d,2)]','EdgeColor','k','LineWidth',1);
%     
end

figs = {fig1, fig2};

for i=[1 2]
    figure(figs{i});
    
    
    axis tight equal off;
    view(1,40);
    camlight;
    view(-2,29);
    
    print(gcf,sprintf('teaser/figure%u.png',i),'-dpng','-r300');
end


fhe = render_mesh.allFaces().halfedge();
for i=1:render_mesh.nv
    render_mesh.v_traits(i).position = render_mesh.v_traits(i).position * scale;
end
v_pos_all = render_mesh.allVertices().position();
fv1 = fhe.from().ind;
fv2 = fhe.to().ind;
fv3 = fhe.next().to().ind;
fv = [fv1;fv2;fv3];
figure('Color','white');
title('Fig. 20 (left)');
hold on;
patch('Vertices',v_pos_all','Faces',fv','EdgeColor','none','FaceColor',light_green,'FaceLighting','gouraud',...
    'SpecularStrength',0.1);
camlight;
axis tight equal off;
view(-2,29);
camlight;

num_planes = size(between_planes,2);
for i=1:num_planes
    a = between_planes(1:3,i);
    b = between_planes(4,i);
    if i==1
        b = b-1.e-5;
    elseif i==num_planes
        b = b+1.e-5;
    end
    [p, he, t, closed] = HETools.traceAllPlaneIntersections(render_mesh, a, b);
    p = p{1};
    plot3(p(1,:),p(2,:),p(3,:),'Color',green,'LineWidth',2);
        
    d2 = [0;0;1];
    d1 = cross(d2, a);

    origin = p(:,1);
    offset = 0.006;
    proj1 = sum((p-origin).*d1,1);
    d1_bounds = [min(proj1)-offset, max(proj1)+offset];
    proj2 = sum((p-origin).*d2,1);
    d2_bounds = [min(proj2), max(proj2)+offset];

    rect = [origin+d1_bounds(1)*d1+d2_bounds(1)*d2,...
        origin+d1_bounds(2)*d1+d2_bounds(1)*d2,...
        origin+d1_bounds(2)*d1+d2_bounds(2)*d2,...
        origin+d1_bounds(1)*d1+d2_bounds(2)*d2];
    patch('Vertices', rect', 'Faces',[1 2 3;1 3 4],'EdgeColor','none','FaceColor',[216 229 216]/255, 'FaceLighting','none','FaceAlpha',0.3);
    patch('Vertices', rect', 'Faces',[1:size(rect,2);[2:size(rect,2), 1]]','EdgeColor',gray,'LineWidth',1);
end
% print(gcf,'teaser/figure3.png','-dpng','-r1200');