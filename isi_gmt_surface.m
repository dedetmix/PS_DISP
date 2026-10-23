function Xf = isi_gmt_surface(X, lon, lat, masker, tension)
% ISI_GMT_SURFACE  Mengisi piksel yang dimasker (mis. ill-conditioned) dengan
% interpolasi GMT surface (continuous curvature splines in tension,
% Smith & Wessel 1990). Membutuhkan GMT 6 terpasang dan bisa dipanggil
% dengan perintah "gmt" dari MATLAB (system).
%
%   X       : npix x ncomp nilai (mis. [dU dE dN])
%   lon,lat : koordinat piksel (grid teratur)
%   masker  : npix x 1 logical, true = piksel yang akan diisi ulang
%   tension : faktor tegangan 0..1 (default 0.35, disarankan GMT untuk data curam)
% Output:
%   Xf      : X dengan piksel masker diganti hasil interpolasi

if nargin < 5 || isempty(tension), tension = 0.35; end
lon = lon(:);  lat = lat(:);  masker = logical(masker(:));
nx = numel(unique(round(lon*1e7)));  ny = numel(unique(round(lat*1e7)));
R  = sprintf('-R%.12f/%.12f/%.12f/%.12f', min(lon), max(lon), min(lat), max(lat));
I  = sprintf('-I%.15g/%.15g', (max(lon)-min(lon))/(nx-1), (max(lat)-min(lat))/(ny-1));  % = grid piksel asli
tmp_in = [tempname '.xyz'];  tmp_out = [tempname '.nc'];

Xf = X;
for q = 1:size(X,2)
    ok = ~masker & ~isnan(X(:,q));
    dlmwrite(tmp_in, [lon(ok) lat(ok) X(ok,q)], 'delimiter', ' ', 'precision', 12);
    cmd = sprintf('gmt surface %s %s %s -T%g -fg -Am -G%s', tmp_in, R, I, tension, tmp_out);
    % MATLAB di Linux memasang LD_LIBRARY_PATH ke pustakanya sendiri (libcurl,
    % libxml2 versi lama) sehingga program luar seperti gmt gagal dimuat.
    % Jalankan gmt dengan LD_LIBRARY_PATH dikosongkan.
    if isunix && ~ismac
        cmd = ['env -u LD_LIBRARY_PATH ' cmd];
    end
    [st, msg] = system(cmd);
    if st ~= 0, error('isi_gmt_surface: GMT gagal:\n%s', msg); end
    vn = {ncinfo(tmp_out).Variables.Name};             % 'lon/lat' (geografis) atau 'x/y'
    if any(strcmp(vn, 'lon')), nxv = 'lon'; nyv = 'lat'; else, nxv = 'x'; nyv = 'y'; end
    gx = double(ncread(tmp_out, nxv));  gy = double(ncread(tmp_out, nyv));
    gz = double(ncread(tmp_out, 'z'))';               % [ny x nx]
    % ambil node grid terdekat (indeks dijepit agar piksel tepi tidak jadi NaN)
    jx = min(max(round((lon(masker) - gx(1)) / (gx(2) - gx(1))) + 1, 1), numel(gx));
    jy = min(max(round((lat(masker) - gy(1)) / (gy(2) - gy(1))) + 1, 1), numel(gy));
    Xf(masker,q) = gz(sub2ind(size(gz), jy, jx));
end
delete(tmp_in);  if exist(tmp_out,'file'), delete(tmp_out); end
end
