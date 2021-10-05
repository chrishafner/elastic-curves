classdef GeometryGenerator < handle
    methods(Static)
        %tapered_end: if true, is tapered to a tip; if numeric positive, is
        %tapered at angle; if negative, width of extrusion is reduced on
        %both sides by that value (so total extrusion width is reduced by
        %double that value)
        %
        function p = generateProfileOutlineLoop(s,half_width,extr_length,tapered_end)
            if nargin < 3
                extr_length = 0;
            end
            if nargin < 4
                tapered_end = false;
            end
            
            if extr_length == 0
                p = [s fliplr(s) s(1); half_width -fliplr(half_width) half_width(1)];
            elseif islogical(tapered_end)
                if ~tapered_end
                    s1 = [s(1)-extr_length, s s(end)+extr_length];
                    hw1 = half_width([1 1:end end]);
                    p = [s1 fliplr(s1) s1(1);
                        hw1 -fliplr(hw1) hw1(1)];
                else
                    s1 = [s(1)-extr_length, s, s(end)+extr_length, fliplr(s), s(1)-extr_length];
                    hw1 = [0, half_width, 0, -fliplr(half_width), 0];
                    p = [s1; hw1];
                end
            elseif isnumeric(tapered_end)
                if tapered_end(1) >= 0
                    if numel(tapered_end) == 1
                        tapered_end = [tapered_end tapered_end];
                    end

                    drop = extr_length * tan(tapered_end);
                    s1 = [s(1)-extr_length, s s(end)+extr_length];
                    hw1 = [half_width(1) - drop(1), half_width, half_width(end) - drop(2)];
                    p = [s1 fliplr(s1) s1(1);
                        hw1 -fliplr(hw1) hw1(1)];
                else
                    if numel(tapered_end) == 1
                        tapered_end = [tapered_end tapered_end];
                    end
                    
                    s1 = [s(1)-extr_length, s([1 1:end end]), s(end)+extr_length];
                    hw1 = [repmat(half_width(1) + tapered_end(1), 1, 2), half_width, repmat(half_width(end) + tapered_end(2), 1, 2)];
                    p = [s1 fliplr(s1) s1(1);
                        hw1 -fliplr(hw1) hw1(1)];
                end
            end
        end
        
        function hw1 = computeSubstripHalfWidth(t1,t2,hw2,K)
            a = (t1+t2)^3-t2^3;
            b = t2^3;
            c_min = b*hw2./K;
            c_max = (a+b)*hw2./K;
            if max(c_min) > min(c_max)
                hw1 = [];
            else
                c = min(c_max);
                hw1 = (1/a) * (K*c - b*hw2);
            end
        end
        
        function F = generateStripFaceList(ncp)
            F = zeros(2*(ncp-1),3);
            F(1:2:end,1) = 1:ncp-1;
            F(1:2:end,2) = ncp+2:2*ncp;
            F(1:2:end,3) = ncp+1:2*ncp-1;
            F(2:2:end,1) = 1:ncp-1;
            F(2:2:end,2) = 2:ncp;
            F(2:2:end,3) = ncp+2:2*ncp;
        end
        
        function loops = generateMicroPerforationProfile(s,width_profile,outer_width,density_width,density_length,pattern_ratio,extrusion_length,smooth_weights)
            if nargin < 8
                smooth_weights = false;
            end
            
            total_length = s(end)-s(1);
            l2 = total_length/(2*density_length*(1+pattern_ratio) + pattern_ratio);
            l1 = pattern_ratio * l2;
            
            num_segs = 4*density_length+1;
            s_break = zeros(1, num_segs+1);
            s_break(1:2:end) = linspace(s(1), s(end)-l1, (num_segs+1)/2);
            s_break(2:2:end) = linspace(s(1)+l1, s(end), (num_segs+1)/2);
            pp_s2ind = PPHelper.makePiecewiseLinear(s,1:length(s));
            s_break_ind = ppval(pp_s2ind, s_break);
            
            if smooth_weights
                wpp_coeff = zeros(num_segs,4);
                wpp_coeff(1:4:end,4) = 1;
                wpp_coeff(2:4:end,1) = 2/l2^3;
                wpp_coeff(2:4:end,2) = -3/l2^2;
                wpp_coeff(2:4:end,4) = 1;
                wpp_coeff(4:4:end,1) = -2/l2^3;
                wpp_coeff(4:4:end,2) = 3/l2^2;
                pp_w = mkpp(s_break, wpp_coeff);
            else
                wpp_val = ones(1,length(s_break));
                wpp_val(3:4:end) = 0;
                wpp_val(4:4:end) = 0;
                pp_w = PPHelper.makePiecewiseLinear(s_break, wpp_val);
            end
            
            h = outer_width - width_profile; % Total hole width at each s
            h1_total = h .* ppval(pp_w, s);
            h1_off = h1_total*0.5/density_width;
            border_thickness = width_profile/(2*density_width);
            h1_first = -outer_width*0.5 + border_thickness + h1_off;
            h1_mid = h1_first .* linspace(1,0,density_width)' + (-h1_first) .* linspace(0,1,density_width)';
            pp_h1_mid = PPHelper.makePiecewiseLinearVector(s,h1_mid);
            
            curves = cell(2*(density_length*(2*density_width-1) + density_width), 1);
            curve_offset = 0;
            
            h1_breaks = [0:4:num_segs-1; 3:4:num_segs+2];
            h1_breaks(1,1) = 1;
            h1_breaks(2,end) = num_segs+1;
            for hi=1:size(h1_breaks,2)
                s_start = s_break(h1_breaks(1,hi));
                h_mid_start = kron(ppval(pp_h1_mid, s_start), [1;1]);
                s_end = s_break(h1_breaks(2,hi));
                h_mid_end = kron(ppval(pp_h1_mid, s_end), [1;1]);
                s_ind_bounds = s_break_ind(h1_breaks(:,hi));
                s_ind_int = ceil(s_ind_bounds(1)):floor(s_ind_bounds(2));
                s_int = s(s_ind_int);
                
                hy = [h_mid_start, ...
                      kron(h1_mid(:,s_ind_int), [1;1]) + repmat(kron(h1_off(s_ind_int), [-1;1]), density_width, 1), ...
                      h_mid_end];
                hx = [s_start, s_int, s_end];
                for ci=1:2*density_width
                    curves{curve_offset+1} = [hx; hy(ci,:)];
                    curve_offset = curve_offset+1;
                end
            end
            
            h2_total = h - h1_total;
            h2_off = h2_total*0.5/(density_width-1);
            h2_mid = 0.5 * (h1_mid(1:end-1,:) + h1_mid(2:end,:));
            pp_h2_mid = PPHelper.makePiecewiseLinearVector(s,h2_mid);
            
            h2_breaks = [2:4:num_segs-3; 5:4:num_segs];
            for hi=1:size(h2_breaks,2)
                s_start = s_break(h2_breaks(1,hi));
                h_mid_start = kron(ppval(pp_h2_mid, s_start), [1;1]);
                s_end = s_break(h2_breaks(2,hi));
                h_mid_end = kron(ppval(pp_h2_mid, s_end), [1;1]);
                
                s_ind_bounds = s_break_ind(h2_breaks(:,hi));
                s_ind_int = ceil(s_ind_bounds(1)):floor(s_ind_bounds(2));
                s_int = s(s_ind_int);
                
                hy = [h_mid_start, ...
                      kron(h2_mid(:,s_ind_int), [1;1]) + repmat(kron(h2_off(s_ind_int), [-1;1]), density_width-1, 1), ...
                      h_mid_end];
                hx = [s_start, s_int, s_end];
                for ci=1:2*(density_width-1)
                    curves{curve_offset+1} = [hx; hy(ci,:)];
                    curve_offset = curve_offset+1;
                end
            end
            
            % concat half-loops
            num_holes = length(curves)/2;
            loops = cell(num_holes+1, 1);
            for i=1:num_holes
                loops{i} = [curves{2*i-1} fliplr(curves{2*i}(:,1:end-1))];
            end
            
            x0 = s(1)-extrusion_length;
            x1 = s(end)+extrusion_length;
            if length(outer_width)==1
                loops{end} = [x0 x1 x1 x0 x0; 0.5*outer_width*[1 1 -1 -1 1]];
            else
                loops{end} = [x0 s x1 x1 fliplr(s) x0 x0; ...
                    0.5*outer_width([1 1:end end]), -0.5*outer_width([end end:-1:1 1]), 0.5*outer_width(1)];
            end
        end
        
        function generatePerforatedProfile(s,width_profile,outer_width,num_switches)
            pp_width = PPHelper.makePiecewiseLinear(s,width_profile);
            hole_profile = outer_width - width_profile;
            pp_s2ind = PPHelper.makePiecewiseLinear(s,1:length(s));
            sbreaks = linspace(s(1), s(end), 4*num_switches+2);
            sbreaks_ind = ppval(pp_s2ind, sbreaks);
            l = sbreaks(2)-sbreaks(1);
            num_segs = 4*num_switches+1;
            ppcoeff = zeros(num_segs,4);
            ppcoeff(1:4:end,4) = 1;
            ppcoeff(2:4:end,1) = 2/l^3;
            ppcoeff(2:4:end,2) = -3/l^2;
            ppcoeff(2:4:end,4) = 1;
            ppcoeff(4:4:end,1) = -2/l^3;
            ppcoeff(4:4:end,2) = 3/l^2;
            pp_weight = mkpp(sbreaks, ppcoeff);
            
            num_h1 = num_switches+1;
            h1_ind = [(0:num_h1-1)*4; (0:num_h1-1)*4+3];
            h1_ind(1,1) = 1;
            h1_ind(2,end) = h1_ind(2,end)-1;
            
            figure;
            hold on;
            hv = Array(10,2);
            for hi=1:num_h1
                
                hv.append([sbreaks(h1_ind(1,hi)), 0]);
                
                s_ind = ceil(sbreaks_ind(h1_ind(1,hi))):floor(sbreaks_ind(h1_ind(2,hi)));
                s_hole = s(s_ind);
                hole_widths = hole_profile(s_ind).*ppval(pp_weight,s_hole);
                hv.append([s_hole', 0.5*hole_widths']);
                hv.append([s_hole', -0.5*hole_widths']);
                
                hv.append([sbreaks(h1_ind(2,hi)), 0]);
                
            end
            
                
            num_h2 = num_switches;
            h2_ind = [(0:num_h2-1)*4+2; (0:num_h2-1)*4+5];
            for hi=1:num_h2
                mp = (2*outer_width - ppval(pp_width, sbreaks(h2_ind(1,hi))))*0.25;
                hv.append([sbreaks(h2_ind(1,hi))*[1;1], [mp;-mp]]);
                
                s_ind = ceil(sbreaks_ind(h2_ind(1,hi))):floor(sbreaks_ind(h2_ind(2,hi)));
                s_hole = s(s_ind);
                weights = ppval(pp_weight,s_hole);
                other_hole_width = hole_profile(s_ind).*weights;
                mid_point = (outer_width + other_hole_width)*0.25;
                hole_width = hole_profile(s_ind) .* (1. - weights) * 0.5;
                
                hv.append([s_hole', mid_point'+0.5*hole_width']);
                hv.append([s_hole', mid_point'-0.5*hole_width']);
                hv.append([s_hole', -mid_point'+0.5*hole_width']);
                hv.append([s_hole', -mid_point'-0.5*hole_width']);
                
                mp = (2*outer_width - ppval(pp_width, sbreaks(h2_ind(2,hi))))*0.25;
                hv.append([sbreaks(h2_ind(2,hi))*[1;1], [mp;-mp]]);
            end
            
            hval = hv.get();
            scatter(hval(:,1), hval(:,2),'.');
            plot(s([1 end end 1 1]), outer_width*0.5*[1 1 -1 -1 1]);
            axis tight equal;
        end

        function [V,N,F] = generateDeformedProfileMesh(p,t,K,extr_length)
            if nargin < 4
                extr_length = 0;
            end
            
            if extr_length ~= 0
                p = [p(:,1)-t(:,1)*extr_length, p, p(:,end)+t(:,end)*extr_length];
                t = t(:,[1 1:end end]);
                K = K(:,[1 1:end end]);
            end

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
    end
end