curves_opt = load('horse/horse_curves.mat');
load('horse/horse_mesh.mat');

curves_opt = curves_opt.curves_opt;
final_scale = 0.005;

% horse_filenames = {'horse/horse_parts_1.obj', 'horse/horse_parts_2b.obj', 'horse/horse_parts_3.obj'};
% for i=1:3
%     meshes(i) = HalfEdgeMesh.importOBJ(horse_filenames{i});
% end

% for mi=1:length(meshes)
%     for i=1:meshes(mi).nv
%         meshes(mi).v_traits(i).position = meshes(mi).v_traits(i).position * final_scale;
%     end
% end
nc = length(curves_opt);
for ci=1:nc
    curves_opt(ci).cp = final_scale * curves_opt(ci).cp;
    curves_opt(ci).gamma = final_scale * curves_opt(ci).gamma;
    curves_opt(ci).origin = final_scale * curves_opt(ci).origin;
end
num_samples = length(curves_opt(1).K);

figure;
% for mi=1:length(meshes)
%     v_pos_all = meshes(mi).allVertices().position();
%     e_all = meshes(mi).allEdges();
%     patch('Vertices', v_pos_all', 'Faces', [e_all.halfedge().from().ind; e_all.halfedge().to().ind]');
% end
% axis equal tight;

v_pos_all = mesh.allVertices().position()*0.5e-2;
e_all = mesh.allEdges();
patch('Vertices', v_pos_all', 'Faces', [e_all.halfedge().from().ind; e_all.halfedge().to().ind]');
axis equal tight;


nc_one_side = nc*0.5;
curve_planes = zeros(4,nc_one_side);
for ci=1:nc_one_side
    curve = curves_opt(ci);
    curve_planes(1:3,ci) = cross(curve.plane(:,1), curve.plane(:,2));
    curve_planes(4,ci) = -dot(curve_planes(1:3,ci), curve.origin);
end

% Compute bisector planes between curve planes, and add planes on the sides
% by reflection
between_planes = zeros(4,nc_one_side+1);
for ci=1:nc_one_side-1
    pl = curve_planes(:,ci)+curve_planes(:,ci+1);
    between_planes(:,ci+1) = pl / norm(pl(1:3));
end
q = curve_planes(:,1);
p1 = between_planes(:,2);
between_planes(:,1) = 2*dot(q(1:3),p1(1:3))*q - p1;
q = curve_planes(:,end);
p2 = between_planes(:,end-1);
between_planes(:,end) = 2*dot(q(1:3),p2(1:3))*q - p2;

fid = fopen('horse/horse_base.txt', 'w');
fs = FsTools(fid);

extr_length = 5.e-3;
thickness = 0.7e-3;

strip_dist = 1.e-3;
lens = zeros(length(curves_opt),1);
hw1_all = zeros(nc, num_samples);
hw2_all = zeros(nc, num_samples);
strip_s = zeros(nc, num_samples);
for ci=1:length(curves_opt)
    curve = curves_opt(ci);
    p = curve.origin + curve.plane * curve.gamma;
    to = p(:,2:end) - p(:,1:end-1);
    seg_lens = sqrt(sum(to.^2,1));
    strip_s(ci,:) = [0 cumsum(seg_lens)];
    lens(ci) = strip_s(ci,end);
    
    normal = cross(curve.plane(:,1), curve.plane(:,2));
    bpi = ci;
    if ci > nc_one_side
        bpi = bpi - nc_one_side;
    end
    
    hw1 = (-between_planes(4,bpi)-sum(p.*between_planes(1:3,bpi),1)) ./ dot(between_planes(1:3,bpi), normal);
    hw2 = (-between_planes(4,bpi+1)-sum(p.*between_planes(1:3,bpi+1),1)) ./ dot(between_planes(1:3,bpi+1), -normal);

    
    
    hw1 = hw1 - strip_dist*0.5;
    hw2 = hw2 - strip_dist*0.5;
    
    hw1_all(ci,:) = hw1;
    hw2_all(ci,:) = hw2;
    
    hw1_ends = hw1([1 end]);
    hw2_ends = hw2([1 end]);
    
    spline = SplineCurve(curve.degree, curve.cp);
    p_ends = p(:,[1 end]);
    tangents = curve.plane * spline.evaluateD([0 spline.t_max]);
    tangents(:,2) = -tangents(:,end);
    tangents = tangents ./ sqrt(sum(tangents.^2,1));
    % Order: tangent, plane normal, thickness
    for ei=[1 2]
        basis = [tangents(:,ei) normal cross(tangents(:,ei), normal)];
        fs.writeCuboid(1e2*[-extr_length -hw2_ends(ei) -0.5*thickness], 1e2*[0 hw1_ends(ei) 0.5*thickness], 1e2*p_ends(:,ei), basis);
    end
    
    ncp = size(p,2);
    v_strip = [p+hw1.*normal, p-hw2.*normal];
    
    f_strip = zeros(2*(ncp-1),3);
    f_strip(1:2:end,1) = 1:ncp-1;
    f_strip(1:2:end,2) = ncp+2:2*ncp;
    f_strip(1:2:end,3) = ncp+1:2*ncp-1;
    f_strip(2:2:end,1) = 1:ncp-1;
    f_strip(2:2:end,2) = 2:ncp;
    f_strip(2:2:end,3) = ncp+2:2*ncp;
    
    patch('Vertices',v_strip', 'Faces', f_strip,'FaceColor','r','EdgeColor','none');
    patch('Vertices', p', 'Faces', [1:size(p,2)-1;2:size(p,2)]', 'EdgeColor','b','LineWidth',2);
end
fclose(fid);
axis tight equal;
title('Rendering of horse model, corresponding to photograph in Fig. 24');
view(20,20);

% Compute width of substrip
t1 = 0.2e-3;
t2 = 0.1e-3;
a = (t1+t2)^3-t2^3;
b = t2^3;
w1_all = zeros(nc, num_samples);
for ci=1:length(curves_opt)
     w2 = hw1_all(ci,:) + hw2_all(ci,:);
     K = curves_opt(ci).K';
     c_min = b*w2./K;
     c_max1 = (a+b)*2*hw1_all(ci,:)./K;
     c_max2 = (a+b)*2*hw2_all(ci,:)./K;
     c_max = min([c_max1 c_max2], [], 1);
     if max(c_min) > min(c_max)
         fprintf('Problem');
     else
         c = min(c_max);
         w1_all(ci,:) = (1/a) * (K*c - b*w2);
     end
end

% Draw strip combos
figure;
for ci=1:length(curves_opt)
    offset = ci*[0;-40.e-3];
    % Thin-paper strip
    v_thin = [[strip_s(ci,:); -hw1_all(ci,:)], fliplr([strip_s(ci,:); hw2_all(ci,:)])] + offset;
    patch('Vertices',v_thin','Faces',[1:size(v_thin,2);[2:size(v_thin,2) 1]]');
    v_thick = [[strip_s(ci,:); -0.5*w1_all(ci,:)], fliplr([strip_s(ci,:); 0.5*w1_all(ci,:)])] + offset;
    patch('Vertices',v_thick','Faces',[1:size(v_thin,2);[2:size(v_thin,2) 1]]','EdgeColor','r');
end
axis tight equal;
title('Fig. 18 (right, 1/2)');

% output strips
extr_length = 5e-3;
hw_reduce = 0.5e-3;
slit_length = 2e-3;

svg_cuts = cell(6*length(curves_opt),1);
for ci=1:length(curves_opt)
    y_offset = ci*(-40.e-3);
    for i=1:2
        x_offset = 0;
        if i==1
            hws = [hw1_all(ci,:); hw2_all(ci,:)];
        else
            hws = repmat(0.5*w1_all(ci,:),2,1);
            x_offset = lens(ci) + 4*extr_length + 5.e-3;
        end
        
        offset = [x_offset;y_offset];
        
        s = strip_s(ci,:);
        s1 = [s(1)-extr_length, s([1 1:end end]), s(end)+extr_length];
        hw1 = [repmat(hws(1,1) - hw_reduce, 1, 2), hws(1,:), repmat(hws(1,end) - hw_reduce, 1, 2)];
        hw2 = [repmat(hws(2,1) - hw_reduce, 1, 2), hws(2,:), repmat(hws(2,end) - hw_reduce, 1, 2)];
        p = [s1 fliplr(s1) s1(1);
            hw1 -fliplr(hw2) hw1(1)];
        svg_cuts{(ci-1)*6+(i-1)*3+1} = p + offset;
        svg_cuts{(ci-1)*6+(i-1)*3+2} = [s(1)-extr_length s(1)-extr_length+slit_length; 0 0] + offset;
        svg_cuts{(ci-1)*6+(i-1)*3+3} = [s(end)+extr_length-slit_length s(end)+extr_length; 0 0] + offset;
    end
end
SvgTools.exportCurves('horse/horse-double-strip.svg', svg_cuts, 1e3/0.352778);

figure;
hold on;
for i=1:length(svg_cuts)
    plot(svg_cuts{i}(1,:), svg_cuts{i}(2,:));
end
axis tight equal;
title('Cutting paths for horse model');
