function [acc_here,pred_pc_val] = fn_decoding_acc_cross(neuron_score, FC_score,max_repetition, cluster_id)




pred_feature = zeros(size(FC_score));
for test_img_idx = 1:1000
    train_img = setdiff(1:1000, test_img_idx);
    train_feature = FC_score(train_img,:);
    train_rsp = neuron_score(:, train_img)';
    test_rsp = neuron_score(:, test_img_idx)';
    [beta, ~, ~] = lscov([train_rsp, ones(size(train_rsp, 1), 1)], train_feature);
    pred_feature(test_img_idx,:)=[test_rsp, ones(size(test_rsp, 1), 1)] * beta;
    fprintf('%d \n', test_img_idx)
end

num_now = max(cluster_id);
num_of_correct = 0;

for ii = 1:num_now
    cluster_pool{ii} = find(cluster_id==ii);
end

for repetition = 1:max_repetition
    to_be_selected = zeros([1,num_now]);
    for cc = 1:num_now
        to_be_selected(cc) = randsample(cluster_pool{cc},1);
    end
    for target_category = 1:num_now
        dist_now = sum((pred_feature(to_be_selected,:)-FC_score(to_be_selected(target_category),:)).^2, 2);
        if(min(dist_now)==dist_now(target_category))
            num_of_correct = num_of_correct + 1;
        end
    end
end

acc_here = num_of_correct./ (num_now*max_repetition);

acc_here = acc_here*100;
