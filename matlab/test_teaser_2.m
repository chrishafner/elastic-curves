load('teaser/teaser_curves.mat');
curves_opt = curves_opt;

scale = 0.03;
num_strips = length(curves_opt);
num_samples = size(curves_opt(1).gamma,2);

for si=1:num_strips
    curves_opt(si).cp = curves_opt(si).cp * scale;
    curves_opt(si).a = curves_opt(si).a * scale;
    curves_opt(si).b = curves_opt(si).b * scale;
    curves_opt(si).gamma = curves_opt(si).gamma * scale;
    curves_opt(si).origin = curves_opt(si).origin * scale;
end

mesh = HalfEdgeMesh.importOBJ('teaser/teaser_sixth.obj');
v_pos_all = mesh.allVertices().position() * scale;
e_all = mesh.allEdges();
% figure;
% patch('Vertices', v_pos_all', 'Faces', [e_all.halfedge().from().ind; e_all.halfedge().to().ind]');
% axis equal tight;
% 
% for si=1:num_strips
%     p = curves_opt(si).origin + curves_opt(si).basis(:,[1 2]) * curves_opt(si).gamma;
%     patch('Vertices',p','Faces',[1:size(p,2)-1;2:size(p,2)]','EdgeColor','r','LineWidth',2);
% end

% Curve planes
curve_planes = zeros(4,num_strips);
for si=1:num_strips
    curve = curves_opt(si);
    normal = curves_opt(si).basis(:,3);
    curve_planes(1:3,si) = normal;
    curve_planes(4,si) = -dot(normal, curve.origin);
end

% Bisector planes and symmetry planes
between_planes = zeros(4,num_strips+1);
between_planes(:,1) = [1;0;0;0];
between_planes(:,end) = [0.5;-sqrt(3)/2;0;0];
for si=1:num_strips-1
    pl = curve_planes(:,si)+curve_planes(:,si+1);
    between_planes(:,si+1) = pl / norm(pl(1:3));
end

strip_dist = 0.5e-3;
lens = zeros(num_strips,1);
hw_all = zeros(num_strips, num_samples);
strip_s = zeros(num_strips, num_samples);
for si=1:num_strips
    curve = curves_opt(si);
    p = curve.origin + curve.basis(:,[1 2]) * curve.gamma;
    to = p(:,2:end) - p(:,1:end-1);
    seg_lens = sqrt(sum(to.^2,1));
    strip_s(si,:) = [0 cumsum(seg_lens)];
    lens(si) = strip_s(si,end);
    
    normal = curve.basis(:,3);
    
    hw1 = (between_planes(4,si)+sum(p.*between_planes(1:3,si),1)) ./ dot(between_planes(1:3,si), normal);
    hw2 = (between_planes(4,si+1)+sum(p.*between_planes(1:3,si+1),1)) ./ dot(between_planes(1:3,si+1), -normal);

    hw1 = hw1 - strip_dist*0.5;
    hw2 = hw2 - strip_dist*0.5;
    hw_all(si,:) = min(hw1, hw2);
end

F = GeometryGenerator.generateStripFaceList(num_samples);
for si=1:num_strips
    curve = curves_opt(si);
    gamma = curve.origin + curve.basis(:,[1 2]) * curve.gamma;
    normal = curve.basis(:,3);
    
    p = [gamma - normal.*hw_all(si,:), gamma + normal.*hw_all(si,:)];
%     patch('Vertices',p','Faces',F,'EdgeColor','none','FaceColor','red');
end

save('teaser/teaser_half_width.mat','hw_all','between_planes');


figure;
strips_cricut = cell(num_strips,1);
extr_length = 5e-3;
for si=1:num_strips
    outer_width = 2*hw_all(si,:);
    K = curves_opt(si).K';
    c_max = min(outer_width ./ K);
    width_profile = K * c_max;
    density_width = 2;
    loops = GeometryGenerator.generateMicroPerforationProfile(strip_s(si,:),width_profile,outer_width,density_width,3,0.2,5e-3);
    offset = [0;-(si-1)*3e-2];
    for li=1:length(loops)
        loops{li} = loops{li} + offset;
        patch('Vertices',loops{li}','Faces',[1:size(loops{li},2)-1;2:size(loops{li},2)]');
    end
    strips_cricut{si} = loops;
end
axis tight equal;
title('Cutting paths for shell model');

SvgTools.exportForCricut('teaser/teaser_cricut.svg', strips_cricut, 1e3/0.352778);


fid = fopen('teaser/teaser_base.txt', 'w');
fs = FsTools(fid);
paper_thickness = 0.7e-3;
hw_tol = 0.5e-3;
for si=1:num_strips
    curve = curves_opt(si);
    spline = SplineCurve(curve.degree, curve.cp);
    p_ends = curve.origin + curve.basis(:,[1 2]) * spline.evaluate([0 spline.t_max]);
    tangents = curve.basis(:,[1 2]) * spline.evaluateD([0 spline.t_max]);
    tangents(:,end) = -tangents(:,end);
    tangents = tangents ./ sqrt(sum(tangents.^2,1));
    normal = curve.basis(:,3);
    hw_ends = hw_all(si,[1 end]);
    % Order: tangent, plane normal, thickness
    for ei=[1 2]
        basis = [tangents(:,ei) normal cross(tangents(:,ei), normal)];
        fs.writeCuboid(1e2*[-extr_length -hw_ends(ei)-hw_tol -0.5*paper_thickness], ...
            1e2*[0 hw_ends(ei)+hw_tol 0.5*paper_thickness], ...
            1e2*p_ends(:,ei), basis);
    end
end
fclose(fid);