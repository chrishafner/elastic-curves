classdef FsTools < handle
    properties
        fid
        var_index = 1
    end
    methods
        function obj = FsTools(fid)
            obj.fid = fid;
        end
        
        function writeCuboid(obj, c1, c2, o, R)
            fprintf(obj.fid, 'fCuboid(context, id + \"%u\", {\"corner1\" : vector([%.7f,%.7f,%.7f])*centimeter, \"corner2\" : vector([%.7f,%.7f,%.7f])*centimeter});\n', ...
                obj.var_index, c1(1), c1(2), c1(3), c2(1), c2(2), c2(3));
            fprintf(obj.fid, 'var v%u = qCreatedBy(id + \"%u\", EntityType.BODY);\n', obj.var_index, obj.var_index);
            fprintf(obj.fid, 'opTransform(context, id + \"%u\", { \"bodies\" : v%u, \"transform\" : {\"linear\" : matrix([[%.7f,%.7f,%.7f],[%.7f,%.7f,%.7f],[%.7f,%.7f,%.7f]]), \"translation\" : vector([%.7f,%.7f,%.7f])*centimeter}});', ...
                obj.var_index+1, obj.var_index, R(1,1), R(1,2), R(1,3), R(2,1), R(2,2), R(2,3), R(3,1), R(3,2), R(3,3), o(1), o(2), o(3));
            obj.var_index = obj.var_index + 2;
        end
    end
end