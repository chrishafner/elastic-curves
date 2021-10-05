classdef HalfEdgeMesh < handle
    properties
        he_from
        he_to
        he_next
        he_prev
        he_twin
        he_e
        he_f
        
        f_he
        e_he
        v_he
        
        v_traits = 0
        he_traits = 0
        
        nv
        ne
        nhe
        nf
    end
    
    methods(Static)
        function obj = importOBJ(filename)
            % Read list of vertices and faces
            fid = fopen(filename);
            lines = textscan(fid,'%s', 'Delimiter', '\n');
            fclose(fid);
            
            lines = lines{1};
            ret = cellfun(@HalfEdgeMesh.classifyObjLine, lines);
            v_lines = lines(ret=='v',:);
            f_lines = lines(ret=='f',:);
            
            % Define vertices -> nv
            v_lines_split = cellfun(@split, v_lines,'UniformOutput',false);
            v = cellfun(@(c) str2double(c(2:end,1)), v_lines_split, 'UniformOutput',false);
            nv = length(v);
            
            % Run a loop of halfedges around every face -> nf, he_from, nhe, f_he, he_f
            f_lines_split = cellfun(@split, f_lines,'UniformOutput',false);
            f = cellfun(@(c) str2double(c(2:end,1)), f_lines_split, 'UniformOutput', false);
            nf = length(f);
            he_from = cell2mat(f);
            nhe = length(he_from);
            f_sizes = cellfun(@length, f);
            f_he = [0; cumsum(f_sizes(1:end-1))]+1;
            he_f = cell2mat(arrayfun(@(i) i*ones(f_sizes(i),1), (1:nf)', 'UniformOutput',false));
            
            % Find he_to
            he_to = [he_from(2:end); 0];
            b = [he_f(2:end)~=he_f(1:end-1); true];
            b_ind = find(b) - f_sizes + 1;
            he_to(b) = he_from(b_ind);
            
            % Find he_prev and he_next
            he_next = [(2:nhe)'; 0];
            he_next(b) = b_ind;
            he_prev = zeros(nhe,1);
            he_prev(he_next) = 1:nhe;
            
            % Find v_he
            v_he = zeros(nv,1);
            v_he(flipud(he_from)) = nhe:-1:1;
            
            % Fine he_twin (0 if there is none)
            twin_lookup = sparse(he_from, he_to, 1:nhe, nv, nv);
            he_twin = full(twin_lookup(sub2ind([nv nv],he_to,he_from)));
            
            % Define edges -> e_he, ne, he_e
            e_chooser = (he_from < he_to | he_twin==0);
            e_he = find(e_chooser);
            ne = length(e_he);
            he_e = (1:nhe)';
            he_e(e_chooser) = 1:ne;
            he_e(~e_chooser) = he_e(he_twin(~e_chooser));
            
            % Store in data structure
            obj = HalfEdgeMesh();
            obj.he_from = he_from';
            obj.he_to = he_to';
            obj.he_next = he_next';
            obj.he_prev = he_prev';
            obj.he_twin = he_twin';
            obj.he_e = he_e';
            obj.he_f = he_f';

            obj.f_he = f_he';
            obj.e_he = e_he';
            obj.v_he = v_he';

            obj.nv = nv;
            obj.ne = ne;
            obj.nhe = nhe;
            obj.nf = nf;
            
            % Modify v_he, so that for boundary vertices, v_he always
            % refers to the boundary he. (For manifold meshes, this
            % halfedge will be unique.)
            he_all = obj.allHalfedges();
            he_bdry = he_all.select(~he_all.hasTwin());
            v_bdry = he_bdry.from();
            obj.v_he(v_bdry.ind) = he_bdry.ind;
            
            % Create vertex traits
            obj.v_traits = struct('position',v);
        end
    end
    
    methods(Static, Access=private)
        function ret = classifyObjLine(str)
            if length(str)>1
                if str(1)=='v' && str(2)==' '
                    ret = 'v';
                elseif str(1)=='f' && str(2)==' '
                    ret = 'f';
                else
                    ret = ' ';
                end
            else
                ret = ' ';
            end
        end
    end
    
    
    methods
        function obj = HalfEdgeMesh()
        end
        
        function ret = allHalfedges(obj)
            ret = HEHalfedge(obj, (1:obj.nhe)');
        end
        
        function ret = allEdges(obj)
            ret = HEEdge(obj, (1:obj.ne)');
        end
        
        function ret = allVertices(obj)
            ret = HEVertex(obj, (1:obj.nv)');
        end
        
        function ret = allFaces(obj)
            ret = HEFace(obj, (1:obj.nf)');
        end
        
        function ret = halfedges(obj, ind)
            ret = HEHalfedge(obj, ind(:));
        end
        
        function ret = edges(obj, ind)
            ret = HEEdge(obj, ind(:));
        end
        
        function ret = vertices(obj, ind)
            ret = HEVertex(obj, ind(:));
        end
        
        function ret = faces(obj, ind)
            ret = HEFace(obj, ind(:));
        end
    end
end