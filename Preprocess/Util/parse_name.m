function [day,subject,gnumber] = parse_name(file_name_here)
dot_loc = find(file_name_here=='_');
day = file_name_here(dot_loc(1)+1:dot_loc(2)-1);
subject_name = file_name_here(dot_loc(2)+1:dot_loc(3)-1);
switch subject_name
    case 'JianJian'
        subject='M1';
    case 'TuTu'
        subject='M5';
    case 'FaCai'
        subject='M2';
    case 'Facai'
        subject='M2';
    case 'ZhuangZhuang'
        subject='M3';
    case 'zhuangzhuang'
        subject='M3';
    case 'MaoDan'
        subject='M4';
    case 'MaoDaN'
        subject='M4';
end
gnumber=file_name_here(end-4);
end