function ok = verify_reference(manifest)
%VERIFY_REFERENCE  Check the preserved reference against its SHA-256 manifest (same digests as `sha256sum -c`).
%   ok = verify_reference                       % reference/baseline_2026-11-20.sha256
%   ok = verify_reference('reference/libuipc_ip_sources.sha256')   % not this: that one lists files in libuipc-ip
% Every line is "<sha256> *<path relative to the manifest's folder>".
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
if nargin < 1, manifest = fullfile(root, 'reference', 'baseline_2026-11-20.sha256'); end
base = fileparts(manifest);
lines = strtrim(splitlines(fileread(manifest)));
lines = lines(~cellfun(@isempty, lines));
bad = {};
for q = 1:numel(lines)
    tok = regexp(lines{q}, '^([0-9a-f]{64}) [ *](.+)$', 'tokens', 'once');
    f = fullfile(base, tok{2});
    if ~exist(f, 'file')
        bad{end+1} = [tok{2} ' (missing)']; %#ok<AGROW>
    elseif ~strcmp(sha256_file(f), tok{1})
        bad{end+1} = [tok{2} ' (digest differs)']; %#ok<AGROW>
    end
end
ok = isempty(bad);
if ok
    fprintf('verify_reference: %d files match %s\n', numel(lines), manifest);
else
    fprintf('verify_reference: %d of %d files FAIL:\n  %s\n', numel(bad), numel(lines), strjoin(bad, '\n  '));
end
end
