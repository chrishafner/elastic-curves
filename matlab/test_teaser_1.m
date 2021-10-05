mesh = HalfEdgeMesh.importOBJ('teaser/teaser_sixth.obj');

v_pos_all = mesh.allVertices().position();
e_all = mesh.allEdges();

% figure;
% patch('Vertices', v_pos_all', 'Faces', [e_all.halfedge().from().ind; e_all.halfedge().to().ind]');
% axis equal tight;

[~,vi_bottom] = min(v_pos_all(2,:));
he1 = mesh.vertices(vi_bottom).halfedge();
he_bdry = HETools.completeBoundaryLoop(he1);
bdry_from = he_bdry.from().position();
he_bdry_start = find(bdry_from(1,:) < 1.e-6, 1, 'last');
he_bdry = mesh.halfedges(he_bdry.ind(he_bdry_start:end));
v_bdry = [he_bdry.from().position(), mesh.halfedges(he_bdry.ind(end)).to().position()];
% patch('Vertices',v_bdry','Faces',[1:size(v_bdry,2)-1;2:size(v_bdry,2)]','EdgeColor','r','LineWidth',2);

num_strips = 10;
v_start = he_bdry.subdividePath(2*num_strips);
v_start = v_start(:,2:2:end);
% patch('Vertices',v_start','Faces',(1:num_strips)','Marker','.','MarkerSize',12,'MarkerEdgeColor','r');


t = linspace(0,1,2*num_strips+1);
t = t(2:2:end);
angle_fun = @(t) (1-t)*(pi/2) + t*(pi/6) - 0.5*t.*(1-t);
angles = angle_fun(t);
normal_angles = angles + pi/2;
plane_a = [cos(normal_angles); sin(normal_angles); zeros(1,num_strips)];
plane_b = -sum(v_start.*plane_a,1);



cyl_p = [0;3.7402];
cyl_r = [5.187;4.808];
cyl1_fun = @(x) 1 - sum(((x([1 2],:) - cyl_p) ./ cyl_r).^2,1);
plane_fun = @(x) 6-x(2,:);
region_fun = @(x) max(cyl1_fun(x), plane_fun(x));

all_curves = cell(num_strips, 1);
all_curves_origin = zeros(3,num_strips);
all_curves_basis = zeros(3,3,num_strips);
all_curves_2d = cell(num_strips, 1);
for si=1:num_strips
    [p, he, t, closed] = HETools.traceAllPlaneIntersections(mesh, plane_a(:,si), plane_b(si));
    p = p{1};
    he = he{1};
    t = t{1};
    closed = closed(1);
    ret_str = HETools.trimPath(region_fun, p, he, t, closed);
    p = ret_str(1).p;
    if ~isempty(ret_str.ep0)
        p = [ret_str.ep0.p, p];
    end
    if ~isempty(ret_str.ep1)
        p = [p, ret_str.ep1.p];
    end
    
    if norm(p(:,1)) > norm(p(:,end))
        p = fliplr(p);
    end
    
%     patch('Vertices',p','Faces',[1:size(p,2)-1;2:size(p,2)]','EdgeColor','r','LineWidth',1);
    
    all_curves{si} = p;
    all_curves_origin(:,si) = p(:,1);
    all_curves_basis(:,:,si) = [[plane_a(2,si); -plane_a(1,si); 0], [0;0;1], -plane_a(:,si)];
    all_curves_2d{si} = all_curves_basis(:,[1 2],si)' * (p - p(:,1));
end


curves_opt = struct('cp',[],'degree',[],'R',[],'K',[],'a',[],'b',[],'t',[],'gamma',[],'origin',[],'basis',[]);
curves_opt(num_strips).cp = [];

num_samples = 2e2;
for si=1:num_strips
    approx = SplineCurveApproximation(all_curves_2d{si},6,4);
%     approx.reg_weight = 1.e-5;
%     approx.max_iter = 250;
    approx.runOptimization(false);
    
    spline = approx.spline;


    SplineCurveTools.optimizeStiffnessRatio(spline, 1.5, 0.1, 0.05, [false;true;false;false], false, 20, 0.05, false, false(2,1));
    t = linspace(0,spline.t_max,num_samples);
    lpopt = LPStiffnessOptimizer(spline, num_samples);
    [K,a,b] = lpopt.optimizeWithInflections();
    fprintf('Strip %u: R = %f\n', si, max(K));
    
    curves_opt(si).cp = spline.cp;
    curves_opt(si).degree = spline.degree;
    curves_opt(si).R = max(K);
    curves_opt(si).K = K;
    curves_opt(si).a = a;
    curves_opt(si).b = b;
    curves_opt(si).t = t;
    curves_opt(si).gamma = spline.evaluate(t);
    curves_opt(si).origin = all_curves_origin(:,si);
    curves_opt(si).basis = all_curves_basis(:,:,si);
end

save('teaser/teaser_curves.mat', 'curves_opt');