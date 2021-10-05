classdef HEHalfedge < matlab.mixin.Copyable
    properties
        mesh
        ind
    end
    
    methods
        function obj = HEHalfedge(mesh, ind)
            obj.mesh = mesh;
            obj.ind = ind;
        end
        
        function ret = twin(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.he_twin(obj.ind));
        end
        
        function ret = next(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.he_next(obj.ind));
        end
        
        function ret = prev(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.he_prev(obj.ind));
        end
        
        function ret = from(obj)
            ret = HEVertex(obj.mesh, obj.mesh.he_from(obj.ind));
        end
        
        function ret = to(obj)
            ret = HEVertex(obj.mesh, obj.mesh.he_to(obj.ind));
        end
        
        function ret = face(obj)
            ret = HEFace(obj.mesh, obj.mesh.he_f(obj.ind));
        end
        
        function ret = edge(obj)
            ret = HEEdge(obj.mesh, obj.mesh.he_e(obj.ind));
        end
        
        function ret = select(obj, sub_ind)
            ret = HEHalfedge(obj.mesh, obj.ind(sub_ind));
        end
        
        function ret = hasTwin(obj)
            ret = obj.twin().ind ~= 0;
        end
        
        % plane <a,p>+b = 0. For each halfedge, gives t, such that
        % <a,(1-t)p+tq>+b = 0. p and q are from and to vertex position
        function ret = intersectPlane(obj, a, b)
            p = obj.from().position();
            q = obj.to().position();
            ret = (sum(a.*p,1)+b) ./ sum(a.*(p-q),1);
        end
        
        % Gives (1-t)p+tq. p and q are from and to vertex position
        function ret = parametricPoint(obj, t)
            ret = (1-t) .* obj.from().position() + t .* obj.to().position();
        end
        
        % Assumes that obj is a directed path of half edges
        % Subdivides the path in num_segs segments of equal length, and
        % and returns the num_segs+1 points bordering them.
        function ret = subdividePath(obj, num_segs)
            p = [obj.from().position, ...
                obj.select(numel(obj.ind)).to().position];
            to = p(:,2:end) - p(:,1:end-1);
            edge_lens = sqrt(sum(to.^2, 1));
            total_len = sum(edge_lens);
            seg_len = total_len / num_segs;
            current_len = 0;
            edge_i = 1;
            edge_t = 0;
            target_len = seg_len;
            ret = zeros(size(p,1), num_segs+1);
            ret(:,1) = p(:,1);
            ret(:,end) = p(:,end);
            output_i = 2;
            while output_i <= num_segs
                edge_rest_length = (1 - edge_t) * edge_lens(edge_i);
                if current_len + edge_rest_length < target_len
                    current_len = current_len + edge_rest_length;
                    edge_i = edge_i + 1;
                    edge_t = 0;
                else
                    del_t = (target_len - current_len) / edge_lens(edge_i);
                    edge_t = edge_t + del_t;
                    ret(:,output_i) = (1-edge_t) * p(:,edge_i) + edge_t * p(:,edge_i+1);
                    current_len = target_len;
                    target_len = target_len + seg_len;
                    output_i = output_i + 1;
                end
            end
        end
    end
end