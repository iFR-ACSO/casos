function copyMfiles(source, destination)
% Copy only m-files.

% get all .m files recursively
files = dir(join([source "**" "*.m"], filesep));

if ~isfolder(destination) 
    % create destination folder 
    mkdir(destination); 
end

for i = 1:length(files)
    % copy each file 
    file_path = fullfile(files(i).folder, files(i).name);
    
    % remove pwd from source path for relative comparison
    relative_path = strrep(files(i).folder, join([pwd string(filesep)],""), "");
    relative_path = strrep(relative_path, source, "");

    if startsWith(relative_path, filesep)
        % remove leading file separator
        relative_path{1}(1) = '';
    end
    
    folder_at_destination = fullfile(destination, relative_path);
    if ~isfolder(folder_at_destination)
        mkdir(folder_at_destination);
    end

    % dest = fullfile(folder_at_destination, files(i).name);
    copyfile(file_path, folder_at_destination);
end

end
