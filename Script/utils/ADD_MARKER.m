function output_img = ADD_MARKER(input_img,x,z,label)

switch label{1}
    case 'Scene'
        color_label = [200,20,200];
    case 'Body'
        color_label = [34,200,0];
    case 'Face'
        color_label = [0,155,248];
    case 'Object'
        color_label = [251,117,0];
    case 'Unknown'
        color_label = [128,128,128];
    case 'Color'
        color_label = [200,158,20];
    otherwise
        color_label = [0,0,0];
end
size = 5;
for xx = -size:1:size
    for zz = -size:1:size
        input_img(x+xx,z+zz,:) = [0,0,0];
    end
end
size = size-1;
for xx = -size:1:size
    for zz = -size:1:size
        input_img(x+xx,z+zz,:) = color_label;
    end
end
output_img = input_img;