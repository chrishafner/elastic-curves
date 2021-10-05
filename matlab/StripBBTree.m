classdef StripBBTree < handle
    % An axis-aligned bounding box tree for a triangle strip, i.e., a list
    % of triangles (1,...,n) such that triangles i and i+1 share an edge
    
    properties
        v
        f
        num_nodes
        nodes
    end
    
    methods
        function obj = StripBBTree(v,f)
            obj.v = v;
            obj.f = f;
            obj.nodes = struct('face_range',[],'bb_min',[],'bb_max',[],'parent',[],'children',[]);
            obj.num_nodes = size(f,2)*2-1;
            obj.nodes(obj.num_nodes).face_range = [];
            
            fx = reshape(obj.v(1,f(:)),3,[]);
            fy = reshape(obj.v(2,f(:)),3,[]);
            fz = reshape(obj.v(3,f(:)),3,[]);
            f_bb_min = [min(fx,[],1); min(fy,[],1); min(fz,[],1)];
            f_bb_max = [max(fx,[],1); max(fy,[],1); max(fz,[],1)];
            
            stack = zeros(4,obj.num_nodes+1);
            stack_ctr = 0;
            stack(:,stack_ctr+1) = [1;nan;1;size(f,2)]; % node_index, parent_index, face_from, face_to
            stack_ctr = stack_ctr + 1;
            next_index = 2;
            while stack_ctr > 0
                val = stack(:,stack_ctr);
                stack_ctr = stack_ctr - 1;
                
                node_index = val(1);
                parent = val(2);
                face_range = val([3 4]);
                obj.nodes(node_index).face_range = face_range;
                obj.nodes(node_index).bb_min = min(f_bb_min(:,face_range(1):face_range(2)), [], 2);
                obj.nodes(node_index).bb_max = max(f_bb_max(:,face_range(1):face_range(2)), [], 2);
                obj.nodes(node_index).parent = parent;
                
                if face_range(2) > face_range(1)
                    mid = floor(0.5*(face_range(1)+face_range(2)));
                    stack(:,stack_ctr+1) = [next_index;node_index;face_range(1);mid];
                    stack(:,stack_ctr+2) = [next_index+1;node_index;mid+1;face_range(2)];
                    obj.nodes(node_index).children = [next_index, next_index+1];
                    next_index = next_index + 2;
                    stack_ctr = stack_ctr + 2;
                end
            end
        end
        
        function [line_segs, t_xi, tree_fi, tree_xi] = intersectTriangle(obj,t)
            bb_min = min(t,[],2);
            bb_max = max(t,[],2);
            
            stack = zeros(obj.num_nodes,1);
            stack_ctr = 0;
            stack(stack_ctr+1) = 1;
            stack_ctr = stack_ctr + 1;
            
            line_segs = Array(10,3);
            t_xi = Array(10,3);
            tree_fi = Array(10,1);
            tree_xi = Array(10,3);
            while stack_ctr > 0
                node = obj.nodes(stack(stack_ctr));
                stack_ctr = stack_ctr - 1;
                if all(node.bb_max > bb_min) && all(bb_max > node.bb_min)
                    if isempty(node.children)
                        f1 = obj.f(:,node.face_range(1));
                        t1 = obj.v(:, f1);
                        [points, xi1, xi2] = StripBBTree.intersectTriangles(t, t1);
                        if ~isempty(points)
                            line_segs.append(points');
                            t_xi.append(xi1');
                            tree_fi.append(node.face_range(1) * ones(size(points,2),1));
                            tree_xi.append(xi2');
                        end
                    else
                        stack(stack_ctr + [1 2]) = node.children;
                        stack_ctr = stack_ctr + 2;
                    end
                end
            end
            line_segs = line_segs.get()';
            t_xi = t_xi.get()';
            tree_fi = tree_fi.get()';
            tree_xi = tree_xi.get()';
        end
    end
    
    methods(Static)
        function [points, xi1, xi2] = intersectTriangles(t1, t2)
            points = zeros(3,2);
            ns = [cross(t1(:,2)-t1(:,1), t1(:,3)-t1(:,1)), cross(t2(:,2)-t2(:,1), t2(:,3)-t2(:,1))];
%             n1 = ns(:,1) / norm(ns(:,1));
%             if abs(dot(n1,t2(:,2)-t2(:,1))) < 1.e-4 && abs(dot(n1,t2(:,3)-t2(:,1))) < 1.e-4
%                 points = zeros(3,0);
%                 xi1 = zeros(3,0);
%                 xi2 = zeros(3,0);
%                 return;
%             end
            
            ind = 0;
            t = cat(3, t1, t2);
            for j=[1 2]
                for i=1:3
                    s = [t(:,1,j)-t(:,3,j), t(:,2,j)-t(:,3,j), t(:,i,3-j)-t(:,mod(i,3)+1,3-j)] \ (t(:,i,3-j)-t(:,3,j));
                    if s(1)>0 && s(1)<1 && s(2)>0 && s(2)<1 && s(1)+s(2)<1 && s(3)>=0 && s(3)<=1
                        ind = ind + 1;
                        points(:,ind) = s(1)*t(:,1,j) + s(2)*t(:,2,j) + (1-s(1)-s(2))*t(:,3,j);
                    end
                end
            end
            points = points(:,1:ind);
            xi = zeros(3,ind,2);
            for i=1:ind
                for j=[1 2]
                    s = [t(:,1,j)-t(:,3,j), t(:,2,j)-t(:,3,j), ns(:,3-j)] \ (points(:,i)-t(:,3,j));
                    xi(:,i,j) = [s(1); s(2); 1-s(1)-s(2)];
                end
            end
            xi1 = xi(:,:,1);
            xi2 = xi(:,:,2);
        end
    end
end