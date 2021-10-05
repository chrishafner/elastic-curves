classdef Array < handle
    properties(GetAccess=public, SetAccess=private)
        len
        capacity
        ncol
    end
    
    properties(Access=private)
        value
    end
    
    methods(Access=public)
        function obj = Array(c,ncol)
            if nargin<1
                c = 10;
            end
            if nargin<2
                ncol = 1;
            end
            obj.ncol = ncol;
            obj.capacity = c;
            obj.value = zeros(obj.capacity, obj.ncol);
            obj.len = 0;
        end
        
        function append(obj, val)
            if (size(val, 1) + obj.len) > obj.capacity
                new_cap = obj.capacity * 2;
                while (size(val, 1) + obj.len) > new_cap
                    new_cap = new_cap*2;
                end
                obj.capacity = new_cap;
                temp = obj.value;
                obj.value = zeros(obj.capacity, obj.ncol);
                obj.value(1:size(temp, 1),:) = temp;
            end
            
            obj.value((obj.len+1):(obj.len+size(val,1)),:) = val;
            obj.len = obj.len + size(val, 1);
        end
        
        function ret = get(obj)
            ret = obj.value(1:obj.len,:);
        end
    end
end