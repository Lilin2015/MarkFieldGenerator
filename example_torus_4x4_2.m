clear
close all
reset_toolkit
warning off

%% parameters
D = load('D_torus.mat').D;   % read the D matrix
K = 2;  % set the alphabet size k

Gen = Generator();  % initialize a generator
Gen.SET_runtime_img = false;    %   display the generation process (only support lattice field)

%% generator
Gen.initial(D,K);
Gen.set_order;
Gen.assign(15000);
            
% save
record = Gen.record;
result = Gen.get_result;