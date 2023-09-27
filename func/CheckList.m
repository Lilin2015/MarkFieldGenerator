classdef CheckList < handle
    
    properties
        data            % checklist, 0-ban, -1-available, 1-selected
        hist_data     % checklist history
        hist_p;         % checklist history pointer
        HIST_LAYER = 5000; % layer limits of checklist history (back to 1 when exceed)
        SET_record_hist = false;
    end
    
    methods
        
        %% construct
        function obj = CheckList(order,K)
            % initial the checklist based on order, the start and end are 0
            obj.data = -ones(length(order),K+2);    % the first two columns are previous and next vertex, the others are labels
            obj.data(order,1)=[0;order(1:end-1)];
            obj.data(order,2)=[order(2:end);0];

            if any(obj.data(:,1)==-1)
                warning('empy index in order'); % some verices are not in order, this is possible if some vertices are preset
            end

            % history
            if obj.SET_record_hist
                obj.hist_data = nan(size(obj.data,1),size(obj.data,2),obj.HIST_LAYER);
                obj.hist_p = 1;
                obj.add_hist;
            end
        end
        
        %% check a value
        function obj=check(obj,index,label)

            assert(~any(obj.data(index,3:end)==1),'some values are already checked');
            obj.data(index,label+3)=1; % 0->3, 1->4, ...

            % clear the ban marks on the right (actually, we only need to clear the next one)
%             v = index;
%             v = obj.next(v);
%             while v~=0
%                 obj.data(v,find([0,0,obj.data(v,3:end)==0]))=-1;
%                 v = obj.next(v);
%             end
            v = obj.next(index);
            if v~=0
                obj.data(v,find([0,0,obj.data(v,3:end)==0]))=-1;
            end
            % record in history
            obj.add_hist;

        end

        %% ban a value
        function obj=ban(obj,index,label)

            assert(obj.data(index,label+3)~=0,'the value is already banned');
            obj.data(index,label+3)=0;
            
            % record in history
            obj.add_hist;

        end
        
        %% read a vertex
        function label = read(obj,index)
            label = find(obj.data(index,3:end)==1)-1;
            assert(length(label)==1,'none or more than one check for a single vertex');
        end

        %% read DOF of a vertex
        function d = DOF(obj,index)
            d = sum(obj.data(index,3:end)==-1,'all');
        end
        
        %% build checklist (in checking order)
        function list = get_checklist(obj)
            v = find(obj.data(:,1)==0,1);
            list = [];
            while true
                list(:,end+1) = [v;obj.data(v,3:end)'];
                v = obj.next(v);
                if v==0
                    return;
                end
            end
        end

        %% get the optimal value
        % S is safety estimation
        function label = opt_label(obj,index,S)
            feasible = (obj.data(index,3:end)==-1);
            
            if ~any(feasible)   % all banned, return empty
                label = [];
                return;
            else
                S(~feasible)=NaN;
                [safety,c] = max(S,[],'omitnan'); % the first one
                
                % comment this if unnecessary (the random one)
%                 c = find(S==max(S,[],'omitnan'));
%                 rr = randperm(length(c));
%                 c = c(rr(1));


                if isinf(safety)    % if the best is zero safety
                    obj.data(index,find(feasible)+2)=0; % ban all
                    label = find(feasible,1,'last')-1;  % select the last one
                    obj.data(index,label+3)=0; % unban the selected one
                else
                    label = c-1;
                end
            end

        end

        %% find the next vertex
        function index = next(obj,index)
            index = obj.data(index,2);
        end
        
        %% find the previous vertex
        function index = previous(obj,index)
            index = obj.data(index,1);
        end

        %% find the unbanned vertex closest to Vi from the left side and in Lconf
        function index = closest(obj,v,Lconf)
            while true
                        if ismember(v,Lconf) && ~all(obj.data(v,3:end)==0)
                            index = v;
                            return;
                        end
                        v = obj.previous(v);
                        if v == 0
                            index = [];
                            return;
                        end
            end
        end

        %% insert Vpre to Vref
        function obj = insert(obj,ref,pre)
            if pre==ref || pre==obj.previous(ref)
                return;
            else
                pre_prev = obj.previous(pre);
                pre_next = obj.next(pre);
                ref_prev = obj.previous(ref);

                % adjust the chain list
                obj.data(pre_prev,2)=pre_next;
                obj.data(pre_next,1)=pre_prev;

                obj.data(ref_prev,2)=pre;
                obj.data(pre,1)=ref_prev;

                obj.data(pre,2)=ref;
                obj.data(ref,1)=pre;
            end
        end
    end

    methods (Access = private)

        %% add history
        function obj = add_hist(obj)
            if obj.SET_record_hist
                obj.hist_data(:,:,obj.hist_p)=obj.data;
                obj.hist_p = obj.hist_p+1;
                if obj.hist_p>obj.HIST_LAYER
                    obj.hist_p=1;
                end
                check_loop(obj);
            end
        end

        %% check if loop exist
        function check_loop(obj)
            loop_step = all(obj.hist_data == obj.hist_data(:,:,obj.hist_p),[1,2]);
            loop_step(obj.hist_p)=0;
            loop_step = find(loop_step, 1);
            if ~isempty(loop_step)
                warning('there is a loop!');
            end
        end
            
    end % methods
end

