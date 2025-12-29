root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));

mkdir(fullfile(root_dir,'Figs','F4','SurfacePLOT_CUT'))
mkdir(fullfile(root_dir,'Figs','F4','SurfacePLOT_CUTVisual'))


input_dir = 'C:\Users\moonl\Desktop\NNN\Figs\F4\SuarfacePLOT';
all_input_img = dir(fullfile(input_dir,'*png'));

[imageData, map, alpha] = imread(fullfile(input_dir,all_input_img(end).name));

for img_idx = 1:length(all_input_img)-1
    file_name_here = all_input_img(img_idx).name;
    img_here = imread(fullfile(input_dir,file_name_here));

    dash_loc = find(file_name_here=='_');
    val =  str2num(file_name_here(dash_loc(end)+1 : end-4));
    if(val<0.2)
        continue
    end
    img_here = im2double(img_here);
    pngimg = im2double(imageData);
    alpha = im2double(alpha);
    img_here = img_here .* (1-alpha) + pngimg.*alpha;

    figure;
    set(gcf,'Position',[5 5 2050 1050])
    imshow(img_here(850:2200,[900:2450,3450:5000],:))
    axis off
    hCbar = colorbar('south',colormap=give_me_orange_bao);
    hCbar.Position(1) = hCbar.Position(1)+0.333;
    hCbar.Position(3) = hCbar.Position(3)/4;
    hCbar.FontSize=15;
    clim([-val,val])
    hCbar.AxisLocation="out";
    hCbar.Ticks = round(linspace(-val,val,5),2);
    hCbar.Label.String='Correlation';
    set_font
    
    figsave(fullfile(root_dir,"Figs/F4/SurfacePLOT_CUT/"), all_input_img(img_idx).name)
    close all

    figure;
    set(gcf,'Position',[5 5 2050 1050])
    imshow(img_here(850:2200,[1750:2450,3450:4150],:))
    axis off
    hCbar = colorbar('south',colormap=give_me_orange_bao);
    hCbar.Position(1) = hCbar.Position(1)+0.18;
    hCbar.Position(3) = hCbar.Position(3)/2;
    hCbar.FontSize=15;
    clim([-val,val])
    hCbar.AxisLocation="out";
    hCbar.Ticks = round(linspace(-val,val,5),2);
    hCbar.Label.String='Correlation';
    set_font
    figsave(fullfile(root_dir,"Figs/F4/SurfacePLOT_CUTVisual/"), all_input_img(img_idx).name)
    close all
end
