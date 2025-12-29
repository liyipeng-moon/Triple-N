function r = give_me_nc(raster, img_idx, intested_img, boot_times)


Image_LOCs = cell(length(intested_img),1);
for ii = 1:length(intested_img)
    Image_LOCs{ii} = find(img_idx==intested_img(ii));
end

d1 = zeros([1,length(intested_img)]);
d2 = zeros([1,length(intested_img)]);
r = zeros([1, boot_times]);
for bb = 1:boot_times
    for ii = 1:length(intested_img)
        LOC_NOW = Image_LOCs{ii};
        data_length = length(LOC_NOW);
        half_points = floor(data_length/2);


        Order_now = randperm(data_length);
        first_half = LOC_NOW(Order_now(1:half_points));
        second_half = LOC_NOW(Order_now(half_points+1:end));

        d1(ii) = mean(raster(first_half));
        d2(ii) = mean(raster(second_half));
    end
    r(bb)=corr(d1',d2');
end
r = mean(r);
r = (2*r)/(1+r);
end