function figsave(result_dir, fig_name)
% Saves the current figure (gcf) in PNG, PDF, and SVG formats to the specified result_dir
% The SVG is saved using a vector print command to ensure high-quality scalable graphics.

set_font
% 
saveas(gcf,fullfile(result_dir, sprintf('%s.png', fig_name)))
saveas(gcf,fullfile(result_dir, sprintf('%s.pdf', fig_name)))
print(gcf,'-vector','-dsvg',fullfile(result_dir, sprintf('%s.svg', fig_name)))