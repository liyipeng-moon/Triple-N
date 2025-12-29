function figsave(result_dir, fig_name)
set_font
% 
saveas(gcf,fullfile(result_dir, sprintf('%s.png', fig_name)))
% saveas(gcf,fullfile(result_dir, sprintf('%s.svg', fig_name)))
saveas(gcf,fullfile(result_dir, sprintf('%s.pdf', fig_name)))
print(gcf,'-vector','-dsvg',fullfile(result_dir, sprintf('%s.svg', fig_name))) % svg