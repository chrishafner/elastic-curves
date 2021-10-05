load('horse/horse_mesh.mat');
load('horse/horse_cross_sections.mat','all_curves');
load('horse/horse_curves.mat','curves_opt');

v_pos_all = mesh.allVertices().position();
e_all = mesh.allEdges();

% figure;
% patch('Vertices', v_pos_all', 'Faces', [e_all.halfedge().from().ind; e_all.halfedge().to().ind]');
% axis equal tight;

fhe = mesh.allFaces().halfedge();
fv1 = fhe.from().ind;
fv2 = fhe.to().ind;
fv3 = fhe.next().to().ind;
fv4 = fhe.prev().from().ind;
fv = [fv1;fv2;fv3;fv4];
figure('Color','white');
patch('Vertices', v_pos_all', 'Faces', fv', 'EdgeColor','none','FaceColor',[0.7 0.7 0.7],...
    'SpecularStrength',0.0,'DiffuseStrength',0.2,'AmbientStrength',0.7);
axis tight equal off;

for i=1:20
    c1 = all_curves(i);
    c2 = all_curves(i+20);
    p1 = c1.origin + c1.plane * c1.p;
    p2 = c2.origin + c2.plane * c2.p;
    p = [p1, p2];
    normal = cross(c1.plane(:,1), c1.plane(:,2));
    d1 = [normal(3);0;-normal(1)];
    d2 = cross(normal,d1);
    
    offset = 2;
    proj1 = sum((p-c1.origin).*d1,1);
    d1_bounds = [min(proj1)-offset, max(proj1)+offset];
    proj2 = sum((p-c1.origin).*d2,1);
    d2_bounds = [min(proj2)-offset, max(proj2)+offset];
    
    rect = [c1.origin+d1_bounds(1)*d1+d2_bounds(1)*d2,...
        c1.origin+d1_bounds(2)*d1+d2_bounds(1)*d2,...
        c1.origin+d1_bounds(2)*d1+d2_bounds(2)*d2,...
        c1.origin+d1_bounds(1)*d1+d2_bounds(2)*d2];
    patch('Vertices', rect', 'Faces',[1 2 3;1 3 4],'EdgeColor','none','FaceColor',[216 229 216]/255, 'FaceLighting','none');
    patch('Vertices', rect', 'Faces',[1:size(rect,2);[2:size(rect,2), 1]]');
    
    patch('Vertices', p1', 'Faces', [1:size(p1,2)-1;2:size(p1,2)]', 'EdgeColor', [52 153 52]/255, 'LineWidth', 2);
    patch('Vertices', p2', 'Faces', [1:size(p2,2)-1;2:size(p2,2)]', 'EdgeColor', [52 153 52]/255, 'LineWidth', 2);
end
view(25,11);
camlight;
title('Fig. 18 (left)');

figure;
hold on;
i_sel = 1:2:17;
for j=1:length(i_sel)
    i = i_sel(j);
    x_offset = mod(j-1, 3)*40;
    y_offset = -floor((j-1)/3)*20;
    
    p = all_curves(i).p;
    p_opt = curves_opt(i).gamma;
    
    plot(p(1,:)+x_offset, p(2,:)+y_offset,'Color',[0.7 0.7 0.7],'LineWidth',2);
    plot(p_opt(1,:)+x_offset, p_opt(2,:)+y_offset,'Color','k','LineWidth',1);
end
axis tight equal off;
title('Fig. 18 (right, 2/2)');