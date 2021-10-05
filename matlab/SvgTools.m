classdef SvgTools < handle
    methods(Static)
        function exportCurves(filename, curves, scale)
            % post-process for output
            p_min = zeros(2,length(curves));
            p_max = p_min;
            for ci = 1:length(curves)
                curves{ci}(2,:) = -curves{ci}(2,:);
                p_min(:,ci) = min(curves{ci},[],2);
                p_max(:,ci) = max(curves{ci},[],2);
            end
            p_min = min(p_min,[],2);
            p_max = max(p_max,[],2);
            sz = p_max - p_min;
            for ci = 1:length(curves)
                curves{ci} = curves{ci} - p_min;
            end
            
            file = fopen(filename, 'w');
            fprintf(file,'<svg width=\"%.4f\" height=\"%.4f\" xmlns=\"http://www.w3.org/2000/svg\">\n', sz(1)*scale, sz(2)*scale);
            for ci = 1:length(curves)
                pstr = sprintf('%.4f,%.4f ', curves{ci}*scale);
                pstr = pstr(1:end-1);
                fprintf(file,'<polyline points=\"');
                fprintf(file,pstr);
                fprintf(file,'\" stroke=\"black\" fill=\"none\"/>\n');
            end
            fprintf(file,'</svg>');
            fclose(file);
        end
        
        % loops is a nested cell array, each cell of loops is one connected
        % object. Each cell of loops{i} is one closed loop of 2d points.
        function exportForCricut(filename, loops, scale)
            f_minmax = @(m) [min(m,[],2) max(m,[],2)];
            res = cell2mat(cellfun(@(c) cell2mat(cellfun(f_minmax, c,'UniformOutput',0)'), loops, 'UniformOutput',0)');
            p_min = min(res(:,1:2:end),[],2);
            p_max = max(res(:,2:2:end),[],2);
            sz = p_max-p_min;
            
            bb_scaled = [p_min(1)*scale, -p_max(2)*scale, sz(1)*scale, sz(2)*scale];
            
            file = fopen(filename, 'w');
            fprintf(file,'<svg version=\"1.1\" xmlns="http://www.w3.org/2000/svg\" x=\"0px\" y=\"0px\" viewBox="%.7f %.7f %.7f %.7f" style="enable-background:new %.7f %.7f %.7f %.7f;" xml:space="preserve">\n', ...
                bb_scaled, bb_scaled);
            
            fprintf(file,'<style type=\"text/css\">\n\t.st0{fill:#00A651;stroke:#000000;stroke-miterlimitstroke-miterlimit:10;}\n</style>\n');
            
            for oi=1:length(loops)
                comp = loops{oi};
                fprintf(file, '\t<path class=\"st0\" d=\"');
                for li=1:length(comp)
                    loop = comp{li} * scale;
                    loop(2,:) = -loop(2,:);
                    fprintf(file, 'M%.7f,%.7f', loop(1,1), loop(2,1));
                    for i=2:size(loop,2)-1
                        fprintf(file, 'L%.7f,%.7f', loop(1,i), loop(2,i));
                    end
                    fprintf(file,'z ');
                end
                fprintf(file, '\"/>\n');
            end
            
            fprintf(file,'</svg>');
            fclose(file);
        end
    end
end