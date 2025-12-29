function orange = give_me_orange_bao()
a = colormap_matplotlib('Oranges');
b = colormap_matplotlib('Blues');
orange = [b([256:-2:1],:);a([1:2:256],:)];
end