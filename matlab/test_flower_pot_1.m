mesh = HalfEdgeMesh.importOBJ('flower-pot/flower-pot.obj');
x_offset = 1.e-3;
for i=1:mesh.nv
    mesh.v_traits(i).position = mesh.v_traits(i).position + x_offset;
end
v = mesh.allVertices().position();

e = mesh.allEdges().halfedge();
e = [e.from().ind; e.to().ind];
% 
% figure(10);
% clf;
% patch('Vertices',v','Faces',e','LineWidth',0.5,'EdgeColor',[0.7 0.7 0.7]);
% hold on;
% axis tight equal vis3d;

he = mesh.allHalfedges();
he_bdry = he.select(~he.hasTwin());
% patch('Vertices',v','Faces',[he_bdry.from().ind; he_bdry.to().ind]',...
%     'LineWidth',2,'EdgeColor','r');
% xlabel('x');
% ylabel('y');

% Extract inner boundary loop
bdry_pos = he_bdry.from().position();
bdry_proj = sum(bdry_pos .* [1;0;-1], 1);
[~,i] = max(bdry_proj);
he_start = he_bdry.select(i);
he_loop = HETools.completeBoundaryLoop(he_start);

num_curves = 16;
angles = linspace(0,2*pi,num_curves+1);
angles = angles(1:end-1);
d1 = [cos(angles); sin(angles); zeros(1,num_curves)];
d2 = repmat([0;0;1], 1, num_curves);
bs = zeros(num_curves, 1);

p_start = zeros(3, num_curves);
curves = cell(num_curves,1);
for ci=1:num_curves
    a = cross(d1(:,ci),d2(:,ci));
%     b = -dot(p_start(:,ci), a);
    b = bs(ci);
    t = he_loop.intersectPlane(a,b);
    b_intersect = t>=0 & t<1;
    t_intersect = t(b_intersect);
    he_intersect = he_loop.select(b_intersect);
    p_intersect = he_intersect.parametricPoint(t_intersect);

    p_proj = sum(p_intersect .* d1(:,ci), 1);
    p_proj(p_proj < 0) = inf;
    [~,i] = min(p_proj);
    he_current = he_intersect.select(i);

    points = Array(10,3);
    t = he_current.intersectPlane(a,b);
    points.append(he_current.parametricPoint(t)');
    while he_current.ind > 0
        t = -1;
        while t <= 0 || t > 1
            he_current = he_current.next();
            t = he_current.intersectPlane(a,b);
        end
        points.append(he_current.parametricPoint(t)');
        he_current = he_current.twin();
    end
    
    curves{ci} = points.get()';
    p_start([1 2], ci) = curves{ci}([1 2],1);
%     patch('Vertices',curves{ci}',...
%         'Faces',[1:size(curves{ci},2)-1;2:size(curves{ci},2)]',...
%         'LineWidth',1);
end

ret = cell(length(curves), 1);
for ci=1:num_curves
    points = curves{ci};
    proj_mat = [d1(:,ci) d2(:,ci)]';
    p2d = proj_mat * (points - p_start(:,ci));
    
    sca = SplineCurveApproximation(p2d, 10, 4);
    sca.runOptimization(false);
    
    ret{ci} = SplineCurveTools.optimizeStiffnessRatio(sca.spline, 2.5, 0.25, 0.05, [true;true;true;true], false);
    
%     ret{ci} = SplineCurveTools.approximateDiscreteCurve(p2d, 4, 6, 2., 0.25, 0.05, [true;true], [true;true], false, false);
    ret{ci}.origin = p_start(:,ci);
    ret{ci}.plane = [d1(:,ci) d2(:,ci)];
    
    p_optim = p_start(:,ci) + [d1(:,ci) d2(:,ci)] * ret{ci}.gamma_sampled;
%     figure(10);
%     patch('Vertices',p_optim','Faces',[1:size(p_optim,2)-1;2:size(p_optim,2)]','EdgeColor',[0.5 0 0.5],'LineWidth',2);
end

clear curves_opt;
curves_opt(num_curves) = ret{num_curves};
for ci=1:num_curves-1
    curves_opt(ci) = ret{ci};
end

% save('flower-pot/curves_R3_ncp10.mat', 'curves_opt');
save('flower-pot/curves_R2.5_ncp10_v2.mat', 'curves_opt');