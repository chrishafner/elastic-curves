classdef ObjTools < handle
    methods(Static)
        function writeOBJ(filename,V,N,F)
            fid = fopen(filename,'w');
            fprintf(fid,'v %.10e %.10e %.10e\r\n',V);
            if isempty(N)
                fprintf(fid,'f %u %u %u\r\n',F);
            else
                fprintf(fid,'vn %.10e %.10e %.10e\r\n',N);
                fprintf(fid,'f %u//%u %u//%u %u//%u\r\n',kron(F,[1;1]));
            end
            fclose(fid);
        end
    end
end