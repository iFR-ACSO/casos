function generate_metadata(path_to_source, version)
% Create CaΣoS meta data.

path_to_meta = join([path_to_source "+casos" "@CasosMeta"], filesep);

%% Create method to query version
% create a method file
lines_version = compose("function v = version(), v = '%s'; end", version);
writelines(lines_version, join([path_to_meta "version.m"], filesep));

%% Create methods to query git commit
% only if inside a git repository
try
    % git repository
    repo = gitrepo;
    last_commit = repo.LastCommit;

    % create method to query git revision
    lines_revision = compose("function v = git_revision(), v = '%s'; end", last_commit.ID);
    writelines(lines_revision, join([path_to_meta "git_revision.m"], filesep));

catch
    % nothing to do
    warning("Not inside a git repository.")
end

end
