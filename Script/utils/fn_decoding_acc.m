function [acc_here,pred_pc_val] = fn_decoding_acc(neuron_score, FC_score,max_repetition,num_series, test_images)


pred_feature = zeros(size(FC_score));
parfor test_img_idx = 1:1000
    train_img = setdiff(1:1000, test_img_idx);
    train_feature = FC_score(train_img,:);
    train_rsp = neuron_score(:, train_img)';
    test_rsp = neuron_score(:, test_img_idx)';
    [beta, ~, ~] = lscov([train_rsp, ones(size(train_rsp, 1), 1)], train_feature);
    pred_feature(test_img_idx,:)=[test_rsp, ones(size(test_rsp, 1), 1)] * beta;
end

FC_score = FC_score(test_images,:);
pred_feature = pred_feature(test_images,:);
for num_now = num_series
    num_of_correct = 0;
    for repetition = 1:max_repetition
        to_be_selected = randperm(length(test_images),num_now);
        dist_now = sum((pred_feature(to_be_selected,:)-FC_score(to_be_selected(1),:)).^2, 2);
        if(min(dist_now)==dist_now(1))
            num_of_correct = num_of_correct + 1;
        end
    end
    acc_here(num_now) = num_of_correct./max_repetition;
end
acc_here = acc_here*100;

pred_pc_val = [];
for PC_here = 1:size(FC_score,2)
    pred_data = pred_feature(:,PC_here);
    actual_data = FC_score(:,PC_here);
    pred_pc_val(PC_here) = corr(actual_data,pred_data);
end