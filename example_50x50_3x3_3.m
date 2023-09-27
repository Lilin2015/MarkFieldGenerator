clear
close all
reset_toolkit
warning off

%% parameters
D = load('D_3_50.mat').D;   % read the D matrix
K = 3;  % set the alphabet size k

Gen = Generator();  % initialize a generator
Gen.SET_runtime_img = true;    %   display the generation process
F_map = zeros(50);

%% generator
color_map = [repmat((0:K-1)'/(K-1),[1,3]);0,1,0;1,0,0;0,0,1]; % different labels are displayed in gray, background is green, conflicts are red and blue


%% generator
Gen.initial(D,K,F_map,color_map);
Gen.set_order;
Gen.assign(15000);
            
% save
record = Gen.record;
result = Gen.get_result;