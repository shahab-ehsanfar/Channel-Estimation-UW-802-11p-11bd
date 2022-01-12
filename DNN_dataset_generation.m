clear variables


fd = 1600;
test_only = 1;
if test_only
   %snr0 = -5:5:30;
   snr0 = -5:5:30;
else
%    snr0 = 35; % do the training for only 35dB
   snr0 = -5:5:30;
end

if test_only
    test_train_ratio = 1;
else
    test_train_ratio = 0.2;
end

for si=1:length(snr0)
    % Load the simulation data
    if length(snr0)==1
        %load(['Raw_DNN_dataset_snr_1']);
        load(['Raw_DNN_dataset_fd_' int2str(fd) '_snr_8']);
    else        
        load(['Raw_DNN_dataset_fd_' int2str(fd) '_snr_' int2str(si)]);
    end
    
    
    if si==1
        all_indices            = randperm(its)';
        testing_data_set_size  = ceil(test_train_ratio*its);
        training_data_set_size = its - testing_data_set_size;
        testing_indices        = all_indices(1:testing_data_set_size);
        training_indices       = all_indices(testing_data_set_size+1:end);
        
        nUSC    = N;
        nSym    = noblk;
        
%         algo    = 'cpUW_noW'; algo_hhat = 'H_hat_out_uwcp_noW_dnn'; algo_h = 'h_genie_UWcp_dnn';
        algo    = 'i11bd_noW'; algo_hhat = 'H_hat_out_11bd_noW_dnn'; algo_h = 'h_genie_CP_11bd_dnn';
        
    end   
    
    Training_DatasetX    = zeros(nUSC, nSym, training_data_set_size);
    Training_DatasetY    = zeros(nUSC, nSym, training_data_set_size);
    Testing_DatasetX     = zeros(nUSC, nSym, testing_data_set_size);
    Testing_DatasetY     = zeros(nUSC, nSym, testing_data_set_size);
    
    Train_X   = zeros(nUSC*2, training_data_set_size * nSym);
    Train_Y   = zeros(nUSC*2, training_data_set_size * nSym);
    Test_X    = zeros(nUSC*2, testing_data_set_size * nSym);
    Test_Y    = zeros(nUSC*2, testing_data_set_size * nSym);
      
    
    for sym_indx = 1:nSym
        % Select the training dataset
        for j = 1:length(training_indices)
            i = training_indices(j);
            H_genie = fft(eval([algo_h '{sym_indx,i}']),N);
            Training_DatasetY(:,sym_indx,j) = H_genie;
            %Training_DatasetX(:,sym_indx,j) = H_hat_out_uwcp_noW_dnn{sym_indx,i};
            Training_DatasetX(:,sym_indx,j) = eval([algo_hhat '{sym_indx,i}']);
        end
        
        % Select the testing data set
        for j = 1:length(testing_indices)
            i = testing_indices(j);
            H_genie = fft(eval([algo_h '{sym_indx,i}']),N);
            Testing_DatasetY(:,sym_indx,j) = H_genie;
%             Testing_DatasetX(:,sym_indx,j) = H_hat_out_uwcp_noW_dnn{sym_indx,i};
            Testing_DatasetX(:,sym_indx,j) = eval([algo_hhat '{sym_indx,i}']);
        end
    end
    
    % Reshape Testing and Training Datasets
    Training_DatasetX_expended = reshape(Training_DatasetX, nUSC, nSym * training_data_set_size);
    Training_DatasetY_expended = reshape(Training_DatasetY, nUSC, nSym * training_data_set_size);
    Testing_DatasetX_expended = reshape(Testing_DatasetX, nUSC, nSym * testing_data_set_size);
    Testing_DatasetY_expended = reshape(Testing_DatasetY, nUSC, nSym * testing_data_set_size);
    
    
    % Complex to Real domain conversion
    Train_X(1:nUSC,:)           = real(Training_DatasetX_expended);
    Train_X(nUSC+1:2*nUSC,:)    = imag(Training_DatasetX_expended);
    Train_Y(1:nUSC,:)           = real(Training_DatasetY_expended);
    Train_Y(nUSC+1:2*nUSC,:)    = imag(Training_DatasetY_expended);
    
    Test_X(1:nUSC,:)              = real(Testing_DatasetX_expended);
    Test_X(nUSC+1:2*nUSC,:)       = imag(Testing_DatasetX_expended);
    Test_Y(1:nUSC,:)              = real(Testing_DatasetY_expended);
    Test_Y(nUSC+1:2*nUSC,:)       = imag(Testing_DatasetY_expended);
    
    % Save training and testing datasets to the DNN_Datasets structure
    DNN_Datasets.('Train_X') =  Train_X;
    DNN_Datasets.('Train_Y') =  Train_Y;
    DNN_Datasets.('Test_X') =  Test_X;
    DNN_Datasets.('Test_Y') =  Test_Y;
    
    DNN_Datasets.('training_indices') = training_indices;
    DNN_Datasets.('testing_indices') = testing_indices;
    DNN_Datasets.('nUSC') = nUSC;
    DNN_Datasets.('nSym') = nSym;
    
    if (length(snr0)==1) 
        si_out = 0; 
    else
        si_out = si;
    end
    
    % Save the DNN_Datasets structure to the specified folder in order to be used later in the Python code
    save(['.\' , algo, '_DNN_Dataset_fd_' int2str(fd) '_snr_' int2str(si_out)],  'DNN_Datasets');
    fprintf('\ndNN dataset are reorganized and stored in \n%s\\%s_DNN_dataset_fd_%s_snr_%s.mat\n',pwd,algo,int2str(fd),int2str(si_out))
    fprintf('Use this file to run the dNN training or testing.\n')
    
end