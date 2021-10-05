classdef HEVertex < matlab.mixin.Copyable
    properties
        mesh
        ind
    end
    
    methods
        function obj = HEVertex(mesh, ind)
            obj.mesh = mesh;
            obj.ind = ind;
        end
        
        function ret = halfedge(obj)
            ret = HEHalfedge(obj.mesh, obj.mesh.v_he(obj.ind));
        end
        
        function ret = select(obj, sub_ind)
            ret = HEVertex(obj.mesh, obj.ind(sub_ind));
        end
        
        function p = position(obj)
            p = [obj.mesh.v_traits(obj.ind).position];
        end
    end
end