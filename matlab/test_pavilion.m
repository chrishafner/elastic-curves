circle_radius = 800.e-3;
arc_angle = pi/3;
ellipse_radii = [75.e-3, 112.5e-3];
num_turns = 1.5;
num_strips = 31;
center_height = 20.e-3;
max_width = 45.e-3;
min_width = 20.e-3;

gamma_func = cell(num_strips,1);
dgamma_func = cell(num_strips,1);
ddgamma_func = cell(num_strips,1);
curv_func = cell(num_strips,1);
t_bounds = zeros(2, num_strips);
strip_frames = zeros(3,3,num_strips);
strip_origins = zeros(3, num_strips);


% Define curves
for si=1:num_strips
    circle_angle = arc_angle*(si-1)/(num_strips-1);
    ellipse_angle = 2*pi*num_turns*(si-1)/(num_strips-1);
    ellipse_basis = [cos(ellipse_angle), -sin(ellipse_angle); ...
        sin(ellipse_angle), cos(ellipse_angle)] .* ellipse_radii;
    
    strip_origins(:,si) = circle_radius * [cos(circle_angle); sin(circle_angle); 0];
    strip_frames(:,:,si) = [[cos(circle_angle); sin(circle_angle); 0], ...
        [0;0;1], ...
        [sin(circle_angle); -cos(circle_angle); 0]];
    
    gamma_func{si} = @(t) ellipse_basis * [cos(t);sin(t)] + [0;center_height];
    dgamma_func{si} = @(t) ellipse_basis * [-sin(t);cos(t)];
    ddgamma_func{si} = @(t) -gamma_func{si}(t);
    curv_func{si} = @(t) curv_func_aux(dgamma_func{si}, ddgamma_func{si}, t);
    
    % Want to find two consecutive intersections with y=0
    t_samp = linspace(0,2*pi,40);
    gamma = gamma_func{si}(t_samp);
    i1 = find(gamma(2,1:end)<=0 & gamma(2,[2:end 1])>0);
    i2 = find(gamma(2,1:end)>0 & gamma(2,[2:end 1])<=0);
    ts = [t_samp(i1); t_samp(i2)];
    if ts(2) < ts(1)
        ts(2) = ts(2) + 2*pi;
    end
    for i=[1 2]
        [ts(i), converged] = NumOpt.newtonsMethod(...
            @(t) sum(gamma_func{si}(t).*[0;1],1), ...
            @(t) sum(dgamma_func{si}(t).*[0;1],1), ...
            ts(i), 1.e-12, 10);
        if ~converged
            fprintf('Root find not converged!');
        end
    end
    
    t_bounds(:,si) = ts;
end

% Compute stiffnesses and widths
use_gravity = true;
special_param = 9.5878;
g = 9.81;
e_gravity = [0; 12*g*special_param];

num_samples = 2e2;
gamma_wp = zeros(3, num_samples, num_strips);
gamma_all = zeros(2, num_samples, num_strips);
half_widths_target = zeros(num_strips, num_samples);
half_widths = zeros(num_strips, num_samples);
K_grav = zeros(num_samples, num_strips);
K_nograv = zeros(num_samples, num_strips);
for si=1:num_strips
    t_samp = linspace(t_bounds(1,si), t_bounds(2,si), num_samples);
    gamma = gamma_func{si}(t_samp);
    gamma_all(:,:,si) = gamma;
    gamma_x = gamma(1,:) + circle_radius;
    half_widths_target(si,:) = gamma_x * tan(pi/num_strips);
    
    gamma_wp(:,:,si) = strip_origins(:,si) + strip_frames(:,[1 2],si)*gamma;
    curv = curv_func{si}(t_samp);
    lpopt = LPStiffnessOptimizer(gamma, curv);
    K_nograv(:,si) = lpopt.optimizeSimple();
    
    to = gamma(:,2:end) - gamma(:,1:end-1);
    seg_lens = sqrt(sum(to.^2,1));
    v_weights = 0.5 * [seg_lens(1) seg_lens(1:end-1)+seg_lens(2:end) seg_lens(end)];
    lpopt.v_weights = v_weights;
    lpopt.e_gravity = e_gravity;
    K_grav(:,si) = lpopt.optimizeWithGravity();
    
    if use_gravity
        K = K_grav(:,si);
    else
        K = K_nograv(:,si);
    end
    half_widths(si,:) = 0.5 * max_width * (K/max(K));
%     half_widths(si,:) = 0.5 * min_width * (K/min(K));
end

% Make plots for paper
paper_plots = true;
if paper_plots
    i_sim = [3 5 8];
    
    figure;
    ax = axes();
    hold on;
    % Simulate with gravity
    for j=1:length(i_sim)
        i = i_sim(j);
        gamma = gamma_all(:,:,i);
        dgamma_ends = dgamma_func{i}(t_bounds(:,i)');
        a01 = atan2(dgamma_ends(2,:), dgamma_ends(1,:));
        [alpha, alpha_start, alpha_end] = AbsoluteAngleElastica.computeAngles(gamma, a01(1), a01(2));
        gamma_to = gamma(:,2:end)-gamma(:,1:end-1);
        seg_lens = sqrt(sum(gamma_to.^2,1));
        elastica = AbsoluteAngleElastica(seg_lens', K_nograv(:,i), alpha_start, alpha_end, gamma(:,1));
        elastica.addEndPointConstraint(gamma(:,end));
        elastica.e_gravity = e_gravity;
        alpha_final = elastica.optimizeWithNewton(alpha);
        elastica.draw(alpha, ax.ColorOrder(j,:));
        elastica.draw(alpha_final, ax.ColorOrder(j,:),'-',2);
    end
    axis tight equal;
    title('Fig. 17 (top left)');
    
    figure;
    ax = axes();
    
    hold on;
    for j=1:length(i_sim)
        i = i_sim(j);
        gamma = gamma_all(:,:,i);
        gamma_to = gamma(:,2:end)-gamma(:,1:end-1);
        seg_lens = sqrt(sum(gamma_to.^2,1));
        s = [0 cumsum(seg_lens)];
        plot(s,K_nograv(:,i)/max(K_nograv(:,i)),'-','LineWidth',1,'Color',ax.ColorOrder(j,:));
        plot(s,K_grav(:,i)/max(K_grav(:,i)),'-','LineWidth',2,'Color',ax.ColorOrder(j,:));
    end
    title('Fig. 17 (top right)');
end

figure;

% build strip geometry
F = GeometryGenerator.generateStripFaceList(num_samples)';
v_strips = zeros(2, num_samples * 2, num_strips);
v_strips_def = zeros(3, num_samples * 2, num_strips);
strip_s = zeros(num_strips, num_samples);
for si=1:num_strips
    normal = strip_frames(:,3,si);
    offsets = half_widths(si,:).*normal;
    v_strips_def(:,:,si) = [gamma_wp(:,:,si) - offsets, gamma_wp(:,:,si) + offsets];
    to = gamma_wp(:,2:end,si) - gamma_wp(:,1:end-1,si);
    seg_lens = sqrt(sum(to.^2,1));
    strip_s(si,:) = [0 cumsum(seg_lens)];
    v_strips(:,:,si) = [[strip_s(si,:);-half_widths(si,:)], [strip_s(si,:);half_widths(si,:)]];
end

% draw frames and curves
indices_draw = 1:num_strips;
for si=indices_draw
    % frame
    frame_scale = 15.e-3;
    for ci=1:3
        col = zeros(3,1);
        col(ci) = 1;
        patch('Vertices', [strip_origins(:,si), strip_origins(:,si) + frame_scale * strip_frames(:,ci,si)]', ...
            'Faces',[1 2],'EdgeColor',col);
    end
    
    % curve
    gamma_world = gamma_wp(:,:,si);
    patch('Vertices', gamma_world', 'Faces', [1:size(gamma_world,2)-1;2:size(gamma_world,2)]', ...
        'EdgeColor','k');
    
    % starting point
    patch('Vertices', gamma_world(:,1)', 'Faces', 1, 'Marker', '.', 'MarkerSize', 6);
    
    % strip
    normal = strip_frames(:,3,si);
    offsets = half_widths(si,:).*normal;
    v = [gamma_world - offsets, gamma_world + offsets];
    
    patch('Vertices', v', 'Faces', F', 'EdgeColor','none','FaceColor','r');
end
axis tight equal vis3d;
camlight;
view(-11,15);
camlight;
view(80,15);
camlight;
title('Rendering of Pavilion, corresponding to photograph in Fig. 1 (left)');

cuts = cell(num_strips,1);
cuts_def = cell(num_strips,1);

% Warnings issued during triangle intersection tests -> intended behavior,
% can ignore them
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

% Compute cuts on every strip in planar space
for si=1:num_strips
    tree = StripBBTree(v_strips_def(:,:,si), F);
    si_others = [si-1 si+1];
    if si==1
        si_others = si+1;
    elseif si==num_strips
        si_others = si-1;
    end
    
    strip_cuts = CellArray();
    strip_cuts_def = CellArray();
    for si_other = si_others
        ls_all = Array(10,3);
        fi_all = Array(10,1);
        xi_all = Array(10,3);
        for fi=1:size(F,2)
            t = v_strips_def(:,F(:,fi),si_other);
            [line_segs, ~, tree_fi, tree_xi] = tree.intersectTriangle(t);
            ls_all.append(line_segs');
            fi_all.append(tree_fi');
            xi_all.append(tree_xi');
        end
        ls_all = ls_all.get()';
        fi_all = fi_all.get()';
        xi_all = xi_all.get()';
        if ~isempty(ls_all)
            inds = connectLineSegmentSoup(ls_all, 1.e-10);
            
            % Get normal pointing towards si_other
            normal = strip_frames(:,3,si);
            if si_other == si+1
                normal = -normal;
            end
            % Orient intersection paths, so they start from boundary
            for ii=1:length(inds)
                p_proj = sum(ls_all(:,inds{ii}) .* normal, 1);
                if p_proj(1) < p_proj(end)
                    inds{ii} = fliplr(inds{ii});
                end
                
                % Extract the first half (in length) of the path
                p = ls_all(:,inds{ii});
                to = p(:,2:end)-p(:,1:end-1);
                seg_lens = sqrt(sum(to.^2,1));
                s = [0 cumsum(seg_lens)];
                total_len = s(end);
                i_last = find(s>total_len/2, 1);
                ss = (total_len/2 - s(i_last-1)) / (s(i_last)-s(i_last-1));
                p_end = (1-ss)*ls_all(:,inds{ii}(i_last-1)) + ss*ls_all(:,inds{ii}(i_last));
                fi_last = -1;
                for fii=[-1 0]
                    fi = fi_all(inds{ii}(i_last + fii));
                    t = v_strips_def(:,F(:,fi),si);
                    n = cross(t(:,2)-t(:,1), t(:,3)-t(:,1));
                    sol = [t(:,1)-t(:,3), t(:,2)-t(:,3), n] \ (p_end-t(:,3));
                    if sol(1)>=0 && sol(1)<=1 && sol(2)>=0 && sol(2)<=1 && sol(1)+sol(2)<=1 && abs(sol(3)) < 1.e-10
                        fi_last = fi;
                        xi_last = [sol(1);sol(2);1-sol(1)-sol(2)];
                        break;
                    end
                end
                
                fi_path = [fi_all(inds{ii}(1:i_last-1)) fi_last];
                xi_path = [xi_all(:,inds{ii}(1:i_last-1)) xi_last];
                v_path_def = v_strips_def(:,F(1,fi_path),si).*xi_path(1,:) + v_strips_def(:,F(2,fi_path),si).*xi_path(2,:) + v_strips_def(:,F(3,fi_path),si).*xi_path(3,:);
                v_path = v_strips(:,F(1,fi_path),si).*xi_path(1,:) + v_strips(:,F(2,fi_path),si).*xi_path(2,:) + v_strips(:,F(3,fi_path),si).*xi_path(3,:);
                
                strip_cuts.append(v_path);
                strip_cuts_def.append(v_path_def);
                
                col = 'b';
                if si_other==si-1
                    col = 'g';
                end
                if ismember(si, indices_draw)
                    patch('Vertices', v_path_def', 'Faces', [1:size(v_path_def,2)-1;2:size(v_path_def,2)]',...
                        'EdgeColor',col,'LineWidth',1);
                    patch('Vertices', v_path_def(:,end)', 'Faces', 1, 'Marker', '.', 'MarkerSize', 6, 'EdgeColor',col);
                end
            end
        end
    end
    cuts{si} = strip_cuts.get();
    cuts_def{si} = strip_cuts_def.get();
end


warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:nearlySingularMatrix');


% For every cut, compute a cut thickness, according to the thickness of the
% material, and the approximate angle of intersection
paper_thickness = 0.2e-3; % 0.2 mm
cut_thicknesses = cell(num_strips,1);
for si=1:num_strips
    xs = Array();
    for ci=1:length(cuts{si})
        cut = cuts{si}{ci};
        if cut(2,1) < cut(2,end)
            si_other = si + 1;
        else
            si_other = si - 1;
        end
        cut_point_def = cuts_def{si}{ci}(:,end);
        [~, i_closest] = min(sum((v_strips_def(:,:,si) - cut_point_def).^2, 1));
        i_closest = mod(i_closest-1, num_samples)+1;
        [~, i_other_closest] = min(sum((v_strips_def(:,:,si_other) - cut_point_def).^2, 1));
        i_other_closest = mod(i_other_closest-1, num_samples)+1;
        
        t = t_bounds(1,si) + (i_closest-1)/(num_samples-1)*(t_bounds(2,si)-t_bounds(1,si));
        t_other = t_bounds(1,si_other) + (i_other_closest-1)/(num_samples-1)*(t_bounds(2,si_other)-t_bounds(1,si_other));
        
        dg = strip_frames(:,[1 2],si) * dgamma_func{si}(t);
        dg = dg / norm(dg);
        dg_other = strip_frames(:,[1 2],si_other) * dgamma_func{si_other}(t_other);
        dg_other = dg_other / norm(dg_other);
        angle = acos(dot(dg,dg_other));
        x = paper_thickness / sin(angle);
        xs.append(x);
    end
    cut_thicknesses{si} = xs.get();
end


% draw planar strips with cuts and export svg
figure;
svg_cuts = CellArray();
extr_length = 10.e-3;
extr_tol = 0; % Width of extrusion is reduced by this much on either side to account for 3d printing inaccuracy
for si=1:num_strips
    offset = -[0;si*55.e-3];
    p = GeometryGenerator.generateProfileOutlineLoop(strip_s(si,:),half_widths(si,:),extr_length, -extr_tol) + offset;
    svg_cuts.append(p(:,[1:end 1]));
    patch('Vertices', p', 'Faces', [1:size(p,2);[2:size(p,2) 1]]');
    for ci=1:length(cuts{si})
        p_cut = cuts{si}{ci}+offset;
        % thicken the cut, and extend slightly
        cut_dir = p_cut(:,end)-p_cut(:,1);
        cut_dir = cut_dir / norm(cut_dir);
        cut_normal = [-cut_dir(2); cut_dir(1)];
        cut_starting_dir = p_cut(:,2) - p_cut(:,1);
        cut_starting_dir = cut_starting_dir / norm(cut_starting_dir);
        p_cut(:,1) = p_cut(:,1) - cut_starting_dir * 2.e-3;
        cut_th = cut_thicknesses{si}(ci) * 0.5;
        cut_path = [p_cut - cut_normal * cut_th, fliplr(p_cut + cut_normal * cut_th)];
        svg_cuts.append(cut_path);
        patch('Vertices',cut_path', 'Faces', [1:size(cut_path,2)-1; 2:size(cut_path,2)]');
    end
end
axis tight equal;
title('Cutting paths for pavilion model (used in Fig. 17, bottom)');
m2pt = 1000/0.3528;
SvgTools.exportCurves('pavilion/strips.svg', svg_cuts.get(), m2pt);

% export cubes for onshape
fid = fopen('pavilion/pavilion_base.txt', 'w');
fs = FsTools(fid);
thickness = 0.7e-3;
width_tol = 0.5e-3; % added to width on either side
for si=1:num_strips
    p_ends = strip_origins(:,si) + strip_frames(:,[1 2],si) * gamma_func{si}(t_bounds(:,si)');
    tangents = strip_frames(:,[1 2],si) * dgamma_func{si}(t_bounds(:,si)');
    tangents(:,2) = -tangents(:,end);
    tangents = tangents ./ sqrt(sum(tangents.^2,1));
    % Order: tangent, plane normal, thickness
    hw_ends = half_widths(si,[1 end]);
    for ei=[1 2]
        basis = [tangents(:,ei) strip_frames(:,3,si) cross(tangents(:,ei), strip_frames(:,3,si))];
        fs.writeCuboid(1e2*[-extr_length -hw_ends(ei)-width_tol -0.5*thickness], ...
            1e2*[0 hw_ends(ei)+width_tol 0.5*thickness], ...
            1e2*p_ends(:,ei), basis);
    end
end
fclose(fid);

function ret = curv_func_aux(dgf, ddgf, t)
    dg = dgf(t);
    ddg = ddgf(t);
    ret = (dg(1,:) .* ddg(2,:) - dg(2,:) .* ddg(1,:)) ./ (sum(dg.^2, 1) .^ 1.5);
end

function ret = connectLineSegmentSoup(ls, tol)
    [~,~,ic] = uniquetol(ls',tol, 'ByRows',true);
    
    ret = CellArray();
    np = size(ls,2);
    while true
        ind_mult = accumarray(ic,1);
        asd = find(ind_mult==1,1);
        if isempty(asd)
            break;
        end
        i_start = find(ic == asd, 1);
        inds = Array();
        inds.append(i_start);
        i = i_start;
        ic(i_start) = np+1;
        while true
            i_next = i+1;
            if mod(i,2)==0
                i_next = i-1;
            end
            inds.append(i_next);
            in_ic = find(ic == ic(i_next));
            ic(in_ic) = np+1;
            if length(in_ic) == 1
                break;
            else
                i = in_ic(in_ic ~= i_next);
            end
        end
        ret.append(inds.get()');
    end
    ret = ret.get();
end