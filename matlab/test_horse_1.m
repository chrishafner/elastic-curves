load('horse/horse_mesh.mat');
v_pos_all = mesh.allVertices().position();
e_all = mesh.allEdges();


fhe = mesh.allFaces().halfedge();
fv1 = fhe.from().ind;
fv2 = fhe.to().ind;
fv3 = fhe.next().to().ind;
fv4 = fhe.prev().from().ind;
fv = [fv1;fv2;fv3;fv4];


final_scale = 0.5;
extr_length = 0.55 / final_scale;


p0 = [-30.42;0;37.3];
q0 = [-26.29;0;64.27];
p1 = [-12.06;0;40.21];
% p2 = [-4.05;0;65.63];
q1 = [-6.2;0;67.09];
p2 = [19.96;0;34.24];
q2 = [14.03;0;65.24];
p3 = [39.57;0;50.45];
q3 = [28.86;0;69.85];
p4 = [49.45;0;62.35];
q4 = [50.18;0;77.99];
p5 = [61.48;0;53.17];
q5 = [66.56;0;60.15];

c1 = [-22.80;0;47.02];
c1r = 10.67;
c1p1 = [-33.42;0;44.51];
c1p2 = [-14.44;0;38.28];

sc1 = [24.55;0;48.72];
sr1 = 15.14;
s1p1 = [20.3;0;34.22];
s1p2 = [39.08;0;45.3];

c2a = [56.54;0;79.18];
c2b = [57.212;0;82.21];
c2r = 2.31;

cp1_normal = cross(c1p2-c1p1,[0;1;0]);
cp1_normal = cp1_normal / norm(cp1_normal);

cyl1_fun = @(x) (c1r^2 - sum((x([1 3],:)-c1([1 3])).^2, 1));
cyl1a_fun = @(x) min(max(cyl1_fun(x), sum(x.*cp1_normal,1)-dot(cp1_normal,c1p1)), sum(x.*[-1;0;0],1)+10);

s1_fun = @(x) (sr1^2 - sum((x-sc1).^2,1));
s1a_fun = @(x) max(max(s1_fun(x), sum(x.*[-1;0;0],1)-dot([-1;0;0],s1p1)), sum(x.*[0;0;1],1)-dot([0;0;1],s1p2));
s1b_fun = @(x) min(s1a_fun(x), sum(x.*[1;0;0],1));

cyl2a_fun = @(x) (sum((x([1 3],:)-c2a([1 3])).^2, 1) - c2r^2);
cyl2b_fun = @(x) (sum((x([1 3],:)-c2b([1 3])).^2, 1) - c2r^2);
cyl2_fun = @(x) min(cyl2a_fun(x), cyl2b_fun(x));

comb_fun = @(x) min(max(cyl1a_fun(x), s1b_fun(x)), cyl2_fun(x));


s = [0;1.3;3;4.2;5.4;6.6];

pp_p = PPHelper.makePiecewiseLinearVector(s,[p0 p1 p2 p3 p4 p5]);
pp_q = PPHelper.makePiecewiseLinearVector(s,[q0 q1 q2 q3 q4 q5]);

ts = linspace(s(1),s(end), 20);

all_curves = struct('p', [], 'origin', [], 'plane', []);
all_curves(2*length(ts)).p = [];

for si=[1 2]
    dist_fun = @(x) min(comb_fun(x), sum(x.*[0;3-si*2;0],1)-extr_length);
    for ti=1:length(ts)
        t = ts(ti);
    %     pt = (1-t)*p1+t*p2;
    %     qt = (1-t)*q1+t*q2;
        pt = ppval(pp_p,t);
        qt = ppval(pp_q,t);
        normal = cross(qt-pt,[0;-1;0]);
        a = normal / norm(normal);
        b = -dot(a,pt);

        [p, he, t, closed] = HETools.traceAllPlaneIntersections(mesh, a, b);
        highest = cellfun(@(q)max(q(3,:)), p);
        [~,i_highest] = max(highest);
        p = p{i_highest};
        he = he{i_highest};
        t = t{i_highest};
        closed = closed(i_highest);
    %     patch('Vertices', p', 'Faces', [1:size(p,2); 2:size(p,2), 1]','EdgeColor','r','LineWidth',2);

        str_trimmed = HETools.trimPath(dist_fun, p, he, t, closed);
        [~,i_trim] = max(arrayfun(@(x) size(x.p,2), str_trimmed));
        str_trimmed = str_trimmed(i_trim);
        p_trim = str_trimmed.p;
        if ~isempty(str_trimmed.ep0)
            p_trim = [str_trimmed.ep0.p p_trim];
        end
        if ~isempty(str_trimmed.ep1)
            p_trim = [p_trim str_trimmed.ep1.p];
        end
        if p_trim(3,1) > p_trim(3,end)
            p_trim = fliplr(p_trim);
        end
        
        e1 = p_trim(:,end) - p_trim(:,1);
        e1 = e1/norm(e1);
        e2 = [0;1;0];
        e2 = e2-dot(e2,e1)*e1;
        e2 = e2/norm(e2);
        plane = [e1,e2];
        origin = p_trim(:,1);
        p2d = plane' * (p_trim - origin);

        ci = (si-1)*length(ts) + ti;
        all_curves(ci).p = p2d;
        all_curves(ci).origin = origin;
        all_curves(ci).plane = plane;
    end
end

save('horse/horse_cross_sections.mat','all_curves');