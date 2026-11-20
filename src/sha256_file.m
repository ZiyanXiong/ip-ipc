function h = sha256_file(f)
%SHA256_FILE  Lower-case hex SHA-256 of a file's bytes (same digest as `sha256sum`).
fid = fopen(f, 'r');
assert(fid > 0, 'sha256_file: cannot open %s', f);
c = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
while true
    b = fread(fid, 2^20, '*uint8');
    if isempty(b), break; end
    md.update(typecast(b, 'int8'));
end
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
