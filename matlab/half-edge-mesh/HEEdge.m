classdef HEEdge < matlab.mixin.Copyable
    properties
        mesh
        ind
    end
    
    methods
        function obj = HEEdge(mesh, ind)
            obj.mesh = mesh;
            obj.ind = ind;
        end
        
        function ret = halfedge(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.e_he(obj.ind));
        end
        
        function ret = select(obj, sub_ind)
            ret = HEEdge(obj.mesh, obj.ind(sub_ind));
        end
    end
end