function metaData = load_meta(meta_file_name)
fprintf('Loading meta data for %s', meta_file_name)
    textData = fileread(meta_file_name);
    lines = strsplit(textData, '\n');
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        
        keyValue = strsplit(line, '=');
        if length(keyValue) ~= 2
            continue;
        end
        
        key = strtrim(keyValue{1});
        value = strtrim(keyValue{2});
        
    
        if any(isstrprop(value, 'digit'))
            if contains(value, '.')
                value = str2double(value);
            else
                value = str2num(value);
            end
        elseif strcmpi(value, 'true')
            value = true;
        elseif strcmpi(value, 'false')
            value = false;
        end
        
        if(key(1)=='~')
            break
        end
        metaData.(key) = value;
    end
end