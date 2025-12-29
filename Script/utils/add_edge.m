function new_img = add_edge(input_img, cc_here, pixsize)
[rows, cols, channels] = size(input_img);
new_rows = rows + 2 * pixsize;
new_cols = cols + 2 * pixsize;

new_img = ones(new_rows, new_cols, channels);
for cc = 1:3
    new_img(:,:,cc) = cc_here(cc);
end
new_img = uint8(new_img);
new_img(pixsize+1:pixsize+rows, pixsize+1:pixsize+cols, :) = input_img;
end