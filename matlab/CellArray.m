classdef CellArray < handle
    properties(GetAccess=public, SetAccess=private)
        len
        capacity
    end
    
    properties(Access=private)
        value
    end
    
    methods(Access=public)
        function obj = CellArray(c)
            if nargin<1
                c = 10;
            end
            obj.capacity = c;
            obj.value = cell(c, 1);
            obj.len = 0;
        end
        
        function append(obj, val)
            if ~iscell(val)
                val = {val};
            end
            if (length(val) + obj.len) > obj.capacity
                new_cap = obj.capacity * 2;
                while (length(val) + obj.len) > new_cap
                    new_cap = new_cap*2;
                end
                obj.capacity = new_cap;
                temp = obj.value;
                obj.value = cell(obj.capacity, 1);
                obj.value(1:length(temp),:) = temp;
            end
            
            obj.value((obj.len+1):(obj.len+length(val)),:) = val;
            obj.len = obj.len + length(val);
        end
        
        function ret = get(obj)
            ret = obj.value(1:obj.len,:);
        end
    end
end