% marker field generator
classdef Generator < handle
    
    properties
        SET_runtime_img = true; % display the generation process
        F_map % a matrix for runtime display, its element num must equal V_num
    end

    properties (GetAccess = public, SetAccess = private)
        D       % the index matrix
        L       % the label matrix
        K       % the alphabet size
        
        R           % a table for search, the i-th row of R contains all the j that D_j (the j-th row of D) contains i, 0 means empty
        Lind       % a table for search, the i-th element of Lind points to the first index of L that corresponds to the i-th label of F
        logU       % a table for fast log calculation 
        
        order   % initial assignment order

        record % assignment history
        
        color_map % label->color, only for runtime display

        V_num   % field size (vertex number of F)
        G_num   % tag number (row number of D)
        G_size   % tag size (the vertex number of G)
    end
    
    methods

        %% construct
        function obj = Generator()
            % empty
        end

        %% initialize
        function obj = initial(obj,D,K,F_map,color_map)

            % take the input parameters
            if nargin < 5
                F_map = [];
                color_map = [];
            end
            obj.assert_D(int64(D));    obj.D = D;
            obj.assert_K(int64(K));     obj.K = K;

            % initialize
            fprintf('\npreparing...');

            % the function did not ask for field (F), thus we can only obtain
            % the size of F by checking the maximum of D. The trivial
            % vertices (not related to any isomorphism) are ignored.
            obj.V_num = max(obj.D,[],'all');
            
            [obj.G_num,obj.G_size] = size(obj.D);

            obj.L = -ones(size(obj.D)); % build the label matrix, having the same size with D

            % take the parameters for display
            if obj.SET_runtime_img
                obj.F_map = F_map;
                obj.color_map = color_map;
                assert(numel(obj.F_map)==obj.V_num,'the element number of F_map must equal V_num, or set SET_runtime_img to false');
                assert(size(obj.color_map,1)==obj.K+3,'the row number of color_map must equal K+3, the last 3 rows are for background, L^conflict, and other conflict, or set SET_runtime_img to false');
            else
                obj.F_map = [];
                obj.color_map = [];
            end

            % reset order, record
            obj.order = [];
            obj.record = [];

            % get R, Lind, logU
            obj.R = zeros(obj.V_num,1);
            for v = 1 : obj.V_num
                [m,~] = find(obj.D==v);
                obj.R(v,1:length(m))=m;
            end
            
            obj.Lind=zeros(obj.V_num,1);
            obj.Lind(flip(obj.D(:)))=flip(1:numel(obj.D));
            
            obj.logU = log(1-(1/obj.K).^(0:obj.G_size));

            fprintf('√');
            
        end
        
        %% precalculate the assignment order
        function set_order(obj)
            % runtime info
            fprintf('\n');
            msg = 0;
            
            L_preserve = obj.L; % preserve the original L, set back latter
            obj.order = zeros(obj.V_num,1);

            % iteratively solve the which problem, and assign the focused
            % vertex by a placeholder
            for i = 1 : obj.V_num
                [measure, f_vertex] = obj.solve_Which;
                if isinf(measure);    break;     end    % if the focused vetex has been assigned, cut off
                
                obj.set(f_vertex,1);          % assign a placeholder
                obj.order(i) = f_vertex;    % record the order
                
                % output info
                msg_c = fprintf([repmat(sprintf('\b'), 1, msg),'setting order...%d/%d'],i,obj.V_num);
                msg = msg_c-msg;
            end
            
            obj.order(obj.order==0)=[]; % cut trivial vertices 

            obj.L = L_preserve; % set back the original L

            fprintf('√');
            
        end

        %% manually set a vertex
        function set(obj,v,label)
            obj.L(obj.D==v) = label;
        end % function

        %% assign label matrix
        function assign(obj,step_limits)

            assert(~isempty(obj.order),'order is empty, set_order first');

            C = CheckList(obj.order,obj.K); % initialize the checklist
            v = obj.order(1);   % get the first focused vertex

            model = 0; % 0 means assignment, 1 means adjustment
            Lconf = []; % the label involved in the current conflict

            % runtime info
            fprintf('\n'); 
            step = 0;
            msg = 0;
            valid_num_max = 0;

            % the loop of assignment and adjustment
            while true
                switch model
                    
                    case 0 % assignment

                        % find the optimal label
                        v_dof = C.DOF(v);
                        if v_dof==1
                            % if v has only one option, the priority does not matter
                            S = zeros(obj.K,1);
                            label = C.opt_label(v,S);
                            obj.set(v,label);   C.check(v,label);
                        elseif v_dof == 0
                            % if v has no option, switch to adjustment without update Lconf
                            model = 1;
                            continue;
                        else
                            S = obj.solve_What(v);
                            label = C.opt_label(v,S);
                            obj.set(v,label);   C.check(v,label);
                            
                            % diplay info
                            step = step + 1;

                            hasL = false(size(obj.L));
                            hasL(all(obj.L~=-1,2),:)=true;
                            valid_num = length(unique(obj.D(hasL)));
                            valid_num_max = max(valid_num_max,valid_num);

                            obj.record(end+1,:)=[step,valid_num_max];
                            msg_c = fprintf([repmat(sprintf('\b'), 1, msg),'assigning...step:%d, valid_num:%d'],step,valid_num_max);
                            msg = msg_c-msg;

                            % break if reach the step limit
                            if step > step_limits
                                break;
                            end
                            
                        end

                        % check conflict
                        [Lconf,Aconf] = obj.find_conflict(v);
                        if isempty(Lconf) % Lconf is empty means no conflict
                            if obj.SET_runtime_img
                                obj.show_Field();
                            end
                            v=C.next(v);
                            if v==0 % v is 0 means no next vertex
                                break;
                            end
                        else % has conflict, switch to adjustment
                            if obj.SET_runtime_img
                                obj.show_Field(Lconf,Aconf);
                            end
                            model = 1;
                            
                            % diplay info
                            obj.record(end+1,:)=[step,valid_num_max];
                            msg_c = fprintf([repmat(sprintf('\b'), 1, msg),'assigning...step:%d, valid_num:%d'],step,valid_num_max);
                            msg = msg_c-msg;

                            continue;
                        end 
                        
                    
                    case 1 % adjustment        
                        vc = C.closest(v,Lconf);   % find the closet involved vertex
                        if isempty(vc)
                            fprintf('\nto root\n');
                            break;
                        end
                        Lconf(Lconf==vc)=[];     % remove the involved vertex before handling it
                        C.insert(v,vc);     % insert the involved vertex to the left of the current vertex
                        v=vc; % the involved vertex becomes the current vertex

                        obj.set(v,-1);   C.ban(v,C.read(v)); % reset the involved vertex, ban its current label
                        model = 0;
                        continue;
                end % switch

            end % while
            
            
            fprintf('√');
            

        end % function

        function results = get_result(obj)
            Lind_m = obj.Lind;
            Lind_m(Lind_m==0)=1;
            results = obj.L(Lind_m);
            results(obj.Lind==0)=-1;
        end

        % show fields
        function show_Field(obj,Lconf,Aconf)
            if isempty(obj.F_map)
                return;
            end
        
            if nargin == 1
                Lconf = [];
                Aconf = [];
            elseif nargin == 2
                Aconf = [];
            end

            % replace vertex labels by colors
            vertex_L = obj.L(obj.Lind);
            vertex_L(vertex_L==-1) = obj.K;
            color = obj.color_map(vertex_L+1,:);
            
            % mask red and blue for Lconf and Aconf vertex
            if ~isempty(Lconf)
                color(Lconf,:)=0.5*color(Lconf,:)+0.5*obj.color_map(obj.K+2,:);
            end
            if ~isempty(Aconf)
                color(Aconf,:)=0.5*color(Aconf,:)+0.5*obj.color_map(obj.K+3,:);
            end
            
            F_show = reshape(color,[size(obj.F_map,1),size(obj.F_map,2),3]);
            imshow(imresize(F_show,[500,500],'nearest'));
            pause(0.001);
        end

    end


    methods (Access = private)

        function assert_D(~,D)  % assert the index matrix
            assert(all(isinteger(D)),'D must be integers');
            assert(all(D>0,'all'),'D must be larger than 0');
            assert(all(~isnan(D),'all'),'D can not be NaN');
            assert(all(~isinf(D),'all'),'D can not be Inf');
        end

        function assert_K(~,K) % assert the alphabet size
            assert(numel(K)==1,'K must be a single value');
            assert(isinteger(K),'K must be integers');
            assert(K>1,'K must be larger than 1');
            assert(~isnan(K),'K can not be NaN');
            assert(~isinf(K),'K can not be Inf');
        end % function

        function [measure, f_vertex] = solve_Which(obj)
            Lexp = [obj.L; nan(1,obj.G_size)];
            Rexp = obj.R;  Rexp(Rexp==0)=obj.G_num+1;
            
            U = sum(Lexp==-1,2);
            U(U==0)=Inf;    % exclude assigned tag
            S = min(U(Rexp),[],2);
                    
            Lind_m = obj.Lind;
            Lind_m(Lind_m==0)=1;
            S([obj.L(Lind_m);0]~=-1)=Inf;
            
            [measure, f_vertex] = min(S);

        end % function

        function S = solve_What(obj, index)

             S = zeros(obj.K,1);    % initialize the safety estimation
             Ri = obj.R(index,:);   % get the indexes of the tags related to Vi
             Ri(Ri==0)=[];            % cut off trivial indexes
            
            for k = 1 : obj.K % calculate the safety of each feasible label
                A = obj.L;
                A(obj.D==index) = k-1; % try the label k

                Srow = []; % the safety of each tag including Vi
        
                Au = (A==-1);
                for r = 1 : length(Ri)
                    ref = A(Ri(r),:);
                    u = Au|(ref==-1); % A or ref is u

                    U = sum(u,2);
                    Spair = obj.logU(U+1);  % U can be 0, then Spair will be -Inf
            
                    Spair(any((~u)&A~=ref,2))=0;
                    Spair(Ri(r))=0;
            
                    Srow(end+1) = sum(Spair);
                end
                S(k) = sum(Srow);
            end
        end % function

        function [Lconf,Aconf] = find_conflict(obj, index)
             Ri = obj.R(index,:);
             Ri(Ri==0)=[];

            % find conflict rows
            isConf = false(obj.G_num,1);
            for i = 1 : length(Ri)
                c_row = Ri(i);
                rep_i = find(all(obj.L==obj.L(c_row,:)&obj.L~=-1,2));
                if length(rep_i)>1
                    isConf(rep_i) = true;
                end
            end

            % find the rows including Vi
            isContain=false(obj.G_num,1);
            isContain(Ri)=1;

            % Lconf is the intersection of both
            isWant = and(isConf,isContain);
            Lconf = unique(obj.D(find(isWant),:));  % conflict vertex in the conflict tags including Vi
            Aconf = unique(obj.D(find(isConf&~isContain),:));  % other conflict vertex
        end % function

    end % methods

end % classdef

