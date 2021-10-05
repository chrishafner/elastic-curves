classdef HEFace < matlab.mixin.Copyable
    properties
        mesh
        ind
    end
    
    methods
        function obj = HEFace(mesh, ind)
            obj.mesh = mesh;
            obj.ind = ind;
        end
        
        function ret = halfedge(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.f_he(obj.ind));
        end

        function ret = select(obj, sub_ind)
            ret = HEFace(obj.mesh, obj.ind(sub_ind));
        end
    end
end