clear
cd C:\Users\moonl\Desktop\NNN
addpath(genpath(pwd));
BrainData = niftiread("NMT_v2.1_sym_SS.nii.gz");
BrainMask = niftiread('NMT_v2.1_sym_brainmask.nii.gz');
NIFTI_DATA = niftiinfo('NMT_v2.1_sym_SS.nii.gz');
BrainData = uint8(BrainData./4);
x0 = size(BrainData,1)/2; % larget is left
y0 = size(BrainData,2)/2; % larger is anterior
z0 = size(BrainData,3)/2; % larger is superior


%%
clear brain_img_to_plot
close all
Area_DATA = readtable('AreaXYZ.xlsx');
area_number = length(Area_DATA.AreaIDX);

AP_series = sort(unique(Area_DATA.A));
img_array = [];

AP_series(1:4)=[];
for AP_IDX = 1:length(AP_series)
    AP_NOW = AP_series(AP_IDX);
    area_here = find(Area_DATA.A==AP_NOW);
    brain_img_here = squeeze(BrainData(:,floor((AP_NOW-NIFTI_DATA.raw.srow_y(end))/0.25),:));
    brain_mask_here = squeeze(BrainMask(:,floor((AP_NOW-NIFTI_DATA.raw.srow_y(end))/0.25),:));

    brain_img_here(find(~brain_mask_here))=255;
    for cc = 1:3
        brain_img_to_plot(:,:,cc) = brain_img_here;
    end

    for aa = 1:length(area_here)
        LR_HERE = Area_DATA.R(area_here(aa));
        LR_coor = floor((LR_HERE-NIFTI_DATA.raw.srow_x(end))/0.25);
        SI_here = Area_DATA.S(area_here(aa));
        SI_coor = floor((SI_here-NIFTI_DATA.raw.srow_z(end))/0.25);
        brain_img_to_plot = ADD_MARKER(brain_img_to_plot,LR_coor,SI_coor,Area_DATA.Label(area_here(aa)));
    end
    for cc = 1:3
        if(AP_NOW<2.5)
            yy_max = 40;
        else
            yy_max = 22;
        end
        for xx = 95:165
            for yy = 19:45
                brain_img_to_plot(xx,yy,cc)=255;
            end
        end
        for xx = 1:size(brain_img_to_plot,1)
            for yy = 1:yy_max
                brain_img_to_plot(xx,yy,cc)=255;
            end
        end
    end
    figure;
    nexttile
    imagesc(imrotate(brain_img_to_plot,90));
    axis equal
    clim([50,255])
    colormap('gray')
    axis off
    hold on
    for aa = 1:length(area_here)
        LR_HERE = Area_DATA.R(area_here(aa));
        LR_coor = floor((LR_HERE-NIFTI_DATA.raw.srow_x(end))/0.25);
        SI_here = Area_DATA.S(area_here(aa));
        SI_coor = floor((SI_here-NIFTI_DATA.raw.srow_z(end))/0.25);
        text(LR_coor-5,200-SI_coor,sprintf('%d',Area_DATA.Subject(area_here(aa))),"FontSize",18,'FontWeight','Bold',FontName='Bold',Color=[1,1,1]);
    end
    fig_here = getframe;
    fig_here = fig_here.cdata(150:305,10:424,:);
    img_array = [img_array,fig_here];
    img_save{AP_IDX}=fig_here;
end
%%
% close all
figure
imshow(img_array)
size_of_one_img = size(fig_here);size_of_one_img = size_of_one_img([1,2]);
img_save{99}=img_save{1}+255;
expected_layout = [reshape(1:24,[6,4])]';
expected_layout(expected_layout>24)=99;
big_img = [];
for yy = 1:size(expected_layout,2)
    row_img = [];
    for xx = 1:size(expected_layout,1)
        if(expected_layout(xx,yy)~=99)
            row_img = [row_img; img_save{expected_layout(xx,yy)}];
        else
            row_img = [row_img; 255*ones(size(img_save{1}))];
        end
    end
    big_img = [big_img, row_img];
end
figure; set(gcf,'Position',[200 200 2200 590])
imshow(big_img)
hold on
for dd = 1:length(AP_series)
    [a,b] = find(expected_layout==dd);
    LOCa_text = (a-1)*size_of_one_img(1);
    LOCb_text = (b-1)*size_of_one_img(2);
    text(170+LOCb_text,125+LOCa_text,sprintf('%.1f mm ',AP_series(dd)),'FontWeight','Bold',"FontSize",12,Color=[0,0,0]);
end
saveas(gcf,'C:\Users\moonl\Desktop\NNN\Figs\F1\MRI.svg')
%%
figure; hold on
set(gcf,'Position',[200 200 700 70])
% set(gca,"Color",[0,0,0])
labels = {'Scene','Body','Face','Object','Color','Unknown'};
for ll = 1:6
    switch labels{ll}
        case 'Scene'
            color_label = [200,20,200];
        case 'Body'
            color_label = [34,200,0];
        case 'Face'
            color_label = [0,155,248];
        case 'Object'
            color_label = [251,117,0];
        case 'Unknown'
            color_label = [128 128 128];
        case 'Color'
            color_label = [200,158,20];
    end
    tt = text(ll,1,labels{ll},'FontWeight','Bold',Color=color_label./255,FontSize=20);
end
axis off
xlim([0.5,7])
ylim([0.8,1.2])
saveas(gcf,'C:\Users\moonl\Desktop\NNN\Figs\F1\LB.svg')