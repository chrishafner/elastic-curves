solver = PlaneCurveSolver();
solver.solve();

fig1 = figure();
clf;
hold on;
plot(solver.parameter_samples, solver.soln.K);
axis([0 2.5 0 5]);
title('Fig. 6 (top right)');

fig2 = figure();
clf;
hold on;
plot(solver.position_samples(1,:), solver.position_samples(2,:));
y = -solver.soln.a/solver.soln.b(2);
plot([-1 1],[y y]);
axis equal;
title('Fig. 6 (top left)');

fig3 = figure();
subplot(3,1,1);
plotProfile(solver.parameter_samples, 0.1*solver.soln.K'/min(solver.soln.K));
axis equal tight off;
sgtitle('Fig. 6 (bottom left)');

[V3,N3,F3] = generateDeformedProfileMesh(solver.position_samples,solver.tangent_samples,0.1*solver.soln.K'/min(solver.soln.K),0.3);
ObjTools.writeOBJ('ellipse-figures/ellipse3.obj',V3,N3,F3);

b = [1 0];
b = b/norm(b);
a = -1;
solver.computeStiffness(a,b);
figure(fig1);
plot(solver.parameter_samples, solver.soln.K/min(solver.soln.K));
figure(fig2);
plot([1 1],[-1 1]);
figure(fig3);
subplot(3,1,2);
plotProfile(solver.parameter_samples, 0.1*solver.soln.K'/min(solver.soln.K));
axis equal tight off;

[V2,N2,F2] = generateDeformedProfileMesh(solver.position_samples,solver.tangent_samples,0.1*solver.soln.K'/min(solver.soln.K),0.3);
ObjTools.writeOBJ('ellipse-figures/ellipse2.obj',V2,N2,F2);

b = [-1 -1];
b = b/norm(b);
a = -0.6;
solver.computeStiffness(a,b);
figure(fig1);
plot(solver.parameter_samples, solver.soln.K/min(solver.soln.K));
figure(fig2);
plot([a*sqrt(2) 0], [0 a*sqrt(2)]);
figure(fig3);
subplot(3,1,3);
plotProfile(solver.parameter_samples, 0.1*solver.soln.K'/min(solver.soln.K));
axis equal tight off;

[V1,N1,F1] = generateDeformedProfileMesh(solver.position_samples,solver.tangent_samples,0.1*solver.soln.K'/min(solver.soln.K),0.3);
ObjTools.writeOBJ('ellipse-figures/ellipse1.obj',V1,N1,F1);


% 3d figure
figure('Color','white');
ax = axes();
hold on;

cols = [126 88 145;
    214 124 83;
    99 159 99]/255;
V = cat(3,V1,V2,V3);
F = cat(3,F1,F2,F3);

for i=1:3
    patch('Vertices',(V(:,:,i) + [0;0.8*(i-1);0])','Faces',F(:,:,i)','EdgeColor','none','FaceColor',cols(i,:),'FaceLighting','gouraud', ...
        'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',0.8);
end
axis tight equal;

xl = ax.XLim;
yl = ax.YLim;
offset = 0.2;
c = [xl(1)-offset, xl(2)+offset;yl(1)-offset, yl(2)+offset;-0.15, 0];
cv = [c(1,1),c(2,1),c(3,1);
    c(1,2),c(2,1),c(3,1);
    c(1,2),c(2,2),c(3,1);
    c(1,1),c(2,2),c(3,1);
    c(1,1),c(2,1),c(3,2);
    c(1,2),c(2,1),c(3,2);
    c(1,2),c(2,2),c(3,2);
    c(1,1),c(2,2),c(3,2)];
cf = [2 7 6; 2 3 7; 3 8 7; 3 4 8; 1 4 3; 1 3 2; 5 7 8; 5 6 7];
patch('Vertices',cv,'Faces',cf,'FaceColor',[1 1 1],'EdgeColor','none');
% scatter3(cv(:,1), cv(:,2), cv(:,3));

for i=1:3
    patch('Vertices',(V(:,:,i).*[1;1;0] + [0;0.8*(i-1);1.e-3])','Faces',F(:,:,i)','EdgeColor','none','FaceColor',[0.6 0.6 0.6],'FaceLighting','none');
end

% patch('Vertices',V1','Faces',F1','EdgeColor','none','FaceColor',cols(1,:),'FaceLighting','gouraud', ...
%     'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',0.8);
% patch('Vertices',(V2 + [0;0.8;0])','Faces',F2','EdgeColor','none','FaceColor',cols(2,:),'FaceLighting','gouraud', ...
%     'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',0.8);
% patch('Vertices',(V3 + [0;1.6;0])','Faces',F3','EdgeColor','none','FaceColor',cols(3,:),'FaceLighting','gouraud', ...
%     'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',0.8);

axis tight equal off;
title('Fig. 6 (bottom right)');



light('Color',[0.5 0.5 0.5],'Position',[0.5 0 1]);
light('Color',[0.5 0.5 0.5],'Position',[0 0.5 1]);
light('Color',[0.2 0.2 0.2],'Position',[0.5 0 -1]);
light('Color',[0.2 0.2 0.2],'Position',[0 0.5 -1]);
view(121,21);

function plotProfile(s,K)
    p = [s fliplr(s) s(1); K -fliplr(K) K(1)];
    plot(p(1,:), p(2,:));
end

function [V,N,F] = generateDeformedProfileMesh(p,t,K,extr_length)
    p = [p(:,1)-t(:,1)*extr_length, p, p(:,end)+t(:,end)*extr_length];
    t = t(:,[1 1:end end]);
    K = K(:,[1 1:end end]);
    
    n = size(p,2);
    normals = [t(2,:); -t(1,:)];
    V = [p(1,:), p(1,:); -K, K; p(2,:), p(2,:)];
    N = [normals(1,:); zeros(1,n); normals(2,:)];
    N = [N, N];
    
    F = zeros(3,2*(n-1));
    F(1,1:2:end) = 1:n-1;
    F(2,1:2:end) = n+1:2*n-1;
    F(3,1:2:end) = n+2:2*n;
    F(1,2:2:end) = 1:n-1;
    F(2,2:2:end) = n+2:2*n;
    F(3,2:2:end) = 2:n;
end