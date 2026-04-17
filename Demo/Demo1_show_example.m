%% This demo code read one session data from session 46
clear
GoodUnit_Dir = 'C:\Users\moonl\Desktop\NNN\Data\Raw\GoodUnit';
Prep_Dir = 'C:\Users\moonl\Desktop\NNN\Data\Processed';
interesred_ses = 46;
load(fullfile(GoodUnit_Dir,'GoodUnit_241112_ZhuangZhuang_NSD1000_LOC_g3.mat'));
load(fullfile(Prep_Dir,'Processed_ses46_241112_M3_3.mat'));
img_idx = meta_data.trial_valid_idx(meta_data.trial_valid_idx~=0); % extract valid trial for this session
psth_t = global_params.PsthRange;
load img_pool.mat
%% how the raster and psth plot of several images
% for the most body selective neuron

% Index for localizer stimuli
face_idx = 1001:1024;
body_idx = 1000+[26:31,43:48,50:61];
obj_idx = setdiff(1025:1072, body_idx);

close all
for example_unit_idx = 1:5
    [xx,neuron_idx] = sort(B_SI,'descend');
    if(isnan(xx(example_unit_idx)))
        continue % some unit show no variance to one category and thus is NAN
    end

    neuron_idx = neuron_idx(example_unit_idx); % unit idx
    response_this_neuron = response_basic(neuron_idx,1:1000); % only extract first 1000
    [~,img_order] = sort(response_this_neuron,'descend');
    example_array = [1:3,500,700,1000]; % select some example index
    col_number = length(example_array);
    figure; set(gcf,'Position',[100 500 1200 500])
    for example_idx = 1:col_number
        example_now = example_array(example_idx); % this order
        img_idx_in_1000 = img_order(example_now); % this img
        trial_wise_location = find(img_idx==img_idx_in_1000); % find trials in Raster data
        raster_here = GoodUnitStrc(neuron_idx).Raster(trial_wise_location,:); % extract raster
        [a,b] = find(raster_here); % plot raster...
        subplot(4,col_number,example_idx)
        scatter(psth_t(b),a,50,'Marker','|','MarkerEdgeColor','k')
        xlim([-30,350])
        ylim([0.5, length(trial_wise_location)+0.5])
        xticks([])
        if(example_idx==1)
            ylabel('#Trial')
        end

        subplot(4,col_number,col_number+example_idx)
        data_to_plot = GoodUnitStrc(neuron_idx).response_matrix_img(img_idx_in_1000,:); % this is avg psth data
        plot(psth_t,data_to_plot,'LineWidth',2,'Color','k')
        xlim([-30,350])
        xlabel('Time (ms)')
        if(example_idx==1)
            ylabel('Firing (Hz)')
            yl = 1.1*max(data_to_plot(:));
        end
        ylim([-1,yl])

        subplot(4,col_number,2*col_number+example_idx)
        imshow(img_pool{img_idx_in_1000}) % show which img
    end

    % show response to localizer stimuli
    subplot(4,col_number,3*col_number+[1:3]); hold on
    data_m = mean(GoodUnitStrc(neuron_idx).response_matrix_img(face_idx,:));
    data_e = std(GoodUnitStrc(neuron_idx).response_matrix_img(face_idx,:))./sqrt(24);
    errorbar(psth_t,data_m,data_e,'LineWidth',1)

    data_m = mean(GoodUnitStrc(neuron_idx).response_matrix_img(body_idx,:));
    data_e = std(GoodUnitStrc(neuron_idx).response_matrix_img(body_idx,:))./sqrt(24);
    errorbar(psth_t,data_m,data_e,'LineWidth',1)

    data_m = mean(GoodUnitStrc(neuron_idx).response_matrix_img(obj_idx,:));
    data_e = std(GoodUnitStrc(neuron_idx).response_matrix_img(obj_idx,:))./sqrt(24);
    errorbar(psth_t,data_m,data_e,'LineWidth',1)

    legend({'Face','Body','Object'},'Box','off','Location','best')
    xlim([-30,350])
    title(sprintf('Unit %d, dprime=%.01f',example_unit_idx,B_SI(neuron_idx)))

    saveas(gcf,sprintf('demo1_unit%d.png',neuron_idx))
end