function res_value = Gm_main(code_result_sp, attributes, n_subspace)
% count the least one between number of discovered attributes and
% meaningful attributes to be the n_subspace

%addpath('H:\Dropbox\manifold_tool\');

%filename = 'min_error_whole_trans_v3_ApAy_32AB';
%filename = 'min_error_whole_trans_v3_fro_ApAy_128AB.mat';

%load(filename);

% PICODES = 1;
% DBC = 2;
% ITQ = 3;
% SPH = 4;
% LSH = 5;
% ATTRIBUTE_GD = 6;

nmethod = size(code_result_sp,1);
 %PROJECTION_KERNEL = @(X,Y) exp(-norm(X * X' - Y * Y','fro')^2/ 1e1);
 GRASMANN_DISTANCE = @Gm_dist_projection;
 metric = GRASMANN_DISTANCE;

%  attributes = attributes_A;
%  code_result_sp = code_result_sp_noise;

 res_table = zeros(nmethod,33);
 
attributes(attributes == 0) = -1;


%[U,D,V] = svd(attributes); %transpose 

%MN_space = U(:,1:n_subspace);

for x_itr = 0
    disp(['iterations : ' num2str(x_itr)]);
    res_table_tmp = zeros(1,nmethod);
    for x_code = 1:nmethod
        DA = code_result_sp{x_code,2 + x_itr};

        DA( DA == 0) = -1;
        
        n_subspace = min(size(DA,1), size(attributes,1));
        
        [U,D,V] = svd(attributes');
        MN_space = U(:,1:n_subspace);
        
        [U,D,V] = svd(DA');
        %[U,D,V] = svd(DA); %transpose 
        DA_space = U(:,1:n_subspace); 

        value = feval(metric,MN_space,DA_space);

        disp([code_result_sp{x_code,1} ' : ' num2str(value)]);
        res_table_tmp(x_code) = -value;
        % Gm_dist(MN_space,DA_space)
    end
    res_table(:,x_itr+1) = res_table_tmp';
    res_value = res_table(:,1);
end