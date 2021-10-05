classdef HETools < handle
    methods(Static)
        function ret = completeBoundaryLoop(he_start)
            if he_start.hasTwin()
                ret = [];
                return;
            end
            
            he_inds = Array(10);
            
            he_current = HEHalfedge(he_start.mesh, he_start.ind);
            while true
                he_inds.append(he_current.ind);
                he_current = he_current.next();
                while he_current.hasTwin()
                    he_current = he_current.twin().next();
                end
                if he_current.ind == he_start.ind
                    break;
                end
            end
            ret = HEHalfedge(he_start.mesh, he_inds.get());
        end
        
        % For every halfedge he, computes the normal of the triangle
        % (he.prev(), he), and stores it as the halfedge trait 'normal'
        function computeHalfedgeNormals(mesh)
            he = mesh.allHalfedges();
            hep = he.prev();
            p0 = hep.from().position();
            p1 = he.from().position();
            p2 = he.to().position();
            n = cross(p1-p0,p2-p0);
            n = n ./ sqrt(sum(n.^2,1));
            if ~isstruct(mesh.he_traits)
                mesh.he_traits = struct('normal', num2cell(n,1)');
            else
                n = num2cell(n,1)';
                [mesh.he_traits(:).normal] = deal(n{:});
            end
        end
        
        % For every halfedge he, computes the angle between the vector
        % (he.prev().twin(), he), and stores it as the halfedge trait 'normal'
        function computeHalfedgeAngles(mesh)
            he = mesh.allHalfedges();
            hep = he.prev();
            p0 = hep.from().position();
            p1 = he.from().position();
            p2 = he.to().position();
            v1 = p0 - p1;
            v1 = v1 ./ sqrt(sum(v1.^2,1));
            v2 = p2 - p1;
            v2 = v2 ./ sqrt(sum(v2.^2,1));
            angles = acos(sum(v1.*v2,1))';
            if ~isstruct(mesh.he_traits)
                mesh.he_traits = struct('angle', num2cell(angles,2));
            else
                a = num2cell(angles,2);
                [mesh.he_traits(:).angle] = deal(a{:});
            end
        end
        
        % Approximates the vertex normal as a weighted sum of halfedge
        % normals, with halfedge angles as weights
        function computeVertexNormals(mesh)
            if ~isstruct(mesh.he_traits)
                HETools.computeHalfedgeNormals(mesh);
                HETools.computeHalfedgeAngles(mesh);
            else
                if ~isfield(mesh.he_traits, 'normal')
                    HETools.computeHalfedgeNormals(mesh);
                end
                if ~isfield(mesh.he_traits, 'angle')
                    HETools.computeHalfedgeAngles(mesh);
                end   
            end
            
            v = mesh.allVertices();
            he_start = v.halfedge();
            v_weights = zeros(1,mesh.nv);
            v_normals = zeros(3,mesh.nv);
            he_current = copy(he_start);
            while ~isempty(he_current.ind)
                v_ind = he_current.from().ind;
                w = [mesh.he_traits(he_current.ind).angle];
                v_weights(v_ind) = v_weights(v_ind) + w;
                n = [mesh.he_traits(he_current.ind).normal];
                v_normals(:,v_ind) = v_normals(:,v_ind) + w .* n;
                
                he_prev = he_current.prev();
                rem = he_prev.hasTwin();
                he_start = he_start.select(rem);
                he_prev = he_prev.select(rem);
                
                he_current = he_prev.twin();
                rem = he_current.ind ~= he_start.ind;
                he_current = he_current.select(rem);
                he_start = he_start.select(rem);
            end
            
            v_normals = v_normals ./ v_weights;
            v_normals = v_normals ./ sqrt(sum(v_normals.^2,1));
            n = num2cell(v_normals,1)';
            [mesh.v_traits(:).normal] = deal(n{:});
        end
        
        % Finds all intersection paths between mesh and the plane
        % <a,x>+b=0. Each connected component is return, along with the
        % information whether it is closed or not. In case of a closed
        % path, the endpoint is not repeated.
        function [p, he, t, closed] = traceAllPlaneIntersections(mesh, a, b)
            he = mesh.allHalfedges();
            t = he.intersectPlane(a,b);
            b_intersect = t>0 & t<1;
            he_intersect = he.select(b_intersect);
            e_intersect = he_intersect.edge();
            e_intersect.ind = unique(e_intersect.ind);
            
            p = CellArray();
            he = CellArray();
            t = CellArray();
            closed = Array();
            
            while ~isempty(e_intersect.ind)
                [p_new, he_new, t_new, closed_new] = HETools.tracePlaneIntersection(a, b, e_intersect.select(1).halfedge(), true);
                p.append(p_new);
                he.append(he_new);
                t.append(t_new);
                closed.append(closed_new);
                
                e_intersect.ind = setdiff(e_intersect.ind, HEHalfedge(mesh, he_new).edge().ind);
            end
            
            p = p.get();
            he = he.get();
            t = t.get();
            closed = closed.get();
        end
        
        % p0 is a cell array of paths. The function returns the index of
        % the path with the vertex closest to point.
        function ind_closest = extractClosestPath(p0, point)
            dist_min = cellfun(@(q) min(sum((q-point).^2, 1)), p0);
            [~, ind_closest] = min(dist_min);
        end
        
        function ret_str = trimPathHalfSpace(a, b, p0, he0, t0, closed0)
            ret_str = HETools.trimPath(@(x)sum(a.*x,1)+b, p0, he0, t0, closed0);
        end
        
        % Intersects the input path the region func(x) >= 0.
        % Returns every connected component of the intersection separately.
        function ret_str = trimPath(func, p0, he0, t0, closed0)
            ep_str = struct('p',[],'he',[],'t',[],'s',[]);
            ret_str = struct('p',[],'he',[],'t',[],'ep0',[],'ep1',[],'closed',[]);
            np = length(t0);
            
            b_in = func(p0) >= 0;
            if all(b_in)
                ret_str.p = p0;
                ret_str.he = he0;
                ret_str.t = t0;
                ret_str.closed = closed0;
            end
            
            i_start = find(~b_in(1:end-1) & b_in(2:end)) + 1;
            if b_in(1)
                i_start = [1 i_start];
            end
            i_end = find(b_in(1:end-1) & ~b_in(2:end));
            if b_in(end)
                i_end = [i_end length(b_in)];
            end
            if closed0 && b_in(1) && b_in(end)
                p0 = [p0,p0];
                he0 = [he0; he0];
                t0 = [t0; t0];
                i_start = i_start(2:end);
                i_end = [i_end(2:end-1), i_end(1) + np];
            end
            
            ret_str(length(i_start)).p = [];
            for i=1:length(i_start)
                i1 = i_start(i);
                i2 = i_end(i);
                ret_str(i).p = p0(:,i1:i2);
                ret_str(i).he = he0(i1:i2);
                ret_str(i).t = t0(i1:i2);
                ret_str(i).closed = false;
                if i_start(i)>1
                    ep_str.s = func(p0(:,i1-1)) / (func(p0(:,i1-1)) - func(p0(:,i1)));
                    ep_str.he = he0([i1-1 i1]);
                    ep_str.t = t0([i1-1 i1]);
                    ep_str.p = (1-ep_str.s)*p0(:,i1-1) + ep_str.s*p0(:,i1);
                    ret_str(i).ep0 = ep_str;
                end
                if i_end(i)<length(t0)
                    ep_str.s = func(p0(:,i2)) / (func(p0(:,i2)) - func(p0(:,i2+1)));
                    ep_str.he = he0([i2 i2+1]);
                    ep_str.t = t0([i2 i2+1]);
                    ep_str.p = (1-ep_str.s)*p0(:,i2) + ep_str.s*p0(:,i2+1);
                    ret_str(i).ep1 = ep_str;
                end
            end
        end
        
        function [p, he, t, is_closed] = tracePlaneIntersection(a, b, he_start, bidirectional)
            if nargin < 4
                bidirectional = false;
            end
            if bidirectional
                [p, he, t, is_closed] = HETools.tracePlaneIntersection(a, b, he_start);
                he_last = HEHalfedge(he_start.mesh, he(end));
                if he_last.hasTwin() || ~he_start.hasTwin()
                    return;
                end
                [p2, he2, t2] = HETools.tracePlaneIntersection(a, b, he_start.twin());
                p = [fliplr(p2(:,2:end)) p];
                he = [flipud(he2(2:end)); he];
                t = [flipud(t2(2:end)); t];
            else
                mesh = he_start.mesh;

                points = Array(10,3);
                he_inds = Array(10,1);
                ts = Array(10,1);
                he_current = HEHalfedge(mesh,he_start.ind);
                t = he_current.intersectPlane(a,b);
                points.append(he_current.parametricPoint(t)');
                he_inds.append(he_current.ind);
                ts.append(t);
                is_closed = false;
                while he_current.ind > 0
                    t = -1;
                    while t < 0 || t > 1
                        he_current = he_current.next();
                        t = he_current.intersectPlane(a,b);
                    end
                    if he_current.twin().ind == he_start.ind
                        is_closed = true;
                        break;
                    end
                    points.append(he_current.parametricPoint(t)');
                    he_inds.append(he_current.ind);
                    ts.append(t);
                    he_current = he_current.twin();
                end

                p = points.get()';
                he = he_inds.get();
                t = ts.get();
            end
        end
    end
end