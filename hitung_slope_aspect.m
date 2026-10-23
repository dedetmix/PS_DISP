function [sE, sN, slope, aspect] = hitung_slope_aspect(dem_file, lon_q, lat_q, sigma, lat0)
% HITUNG_SLOPE_ASPECT  Gradien lereng, slope dan aspect dari DEM (GMT .grd / netCDF)
% pada titik-titik (lon_q, lat_q). Algoritma sama dengan simulasi Python
% (gaussian smoothing sigma piksel, lalu gradien beda-tengah).
%
%   sE     : dz/dE (tanpa satuan)  -> koefisien SPF untuk dE
%   sN     : dz/dN (tanpa satuan)  -> koefisien SPF untuk dN
%   slope  : kemiringan lereng (derajat)
%   aspect : azimut arah TURUN lereng dari utara, searah jarum jam (derajat)
%
% dem_file : grid geografis (bujur/lintang dalam derajat, elevasi dalam meter);
%            nama variabel lon/lat/z atau x/y/z (keduanya dikenali)
% sigma    : smoothing Gaussian dalam piksel DEM (default 2)
% lat0     : lintang acuan konversi derajat bujur -> meter
%            (default: rata-rata lintang titik yang diminta)

lon_q = lon_q(:);  lat_q = lat_q(:);
if nargin < 4 || isempty(sigma), sigma = 2; end
if nargin < 5 || isempty(lat0),  lat0  = mean(lat_q); end

% --- cek file & nama variabel ---
if ~exist(dem_file, 'file')
    error(['hitung_slope_aspect: file DEM "%s" tidak ditemukan di folder kerja (%s).\n' ...
           'Salin DEM area studi ke folder ini, isi dem_file dengan path lengkap, ' ...
           'atau pakai sumber_lereng = ''data''.'], dem_file, pwd);
end
info = ncinfo(dem_file);
vn = {info.Variables.Name};
nx_ = cariNama(vn, {'lon','x','longitude'});
ny_ = cariNama(vn, {'lat','y','latitude'});
nz_ = cariNama(vn, {'z','elevation','Band1','height'});
if isempty(nx_) || isempty(ny_) || isempty(nz_)
    error('hitung_slope_aspect: variabel lon/lat/z tidak dikenali di %s (ada: %s).', ...
          dem_file, strjoin(vn, ', '));
end

x = double(ncread(dem_file, nx_));  x = x(:).';
y = double(ncread(dem_file, ny_));  y = y(:);
if min(lon_q) < min(x) || max(lon_q) > max(x) || min(lat_q) < min(y) || max(lat_q) > max(y)
    error(['hitung_slope_aspect: DEM tidak mencakup area data.\n' ...
           '  DEM : bujur %.4f..%.4f, lintang %.4f..%.4f\n' ...
           '  data: bujur %.4f..%.4f, lintang %.4f..%.4f\n' ...
           'Pastikan DEM milik area studi yang sama dan dalam derajat geografis (bukan UTM).'], ...
          min(x), max(x), min(y), max(y), min(lon_q), max(lon_q), min(lat_q), max(lat_q));
end
flip_x = x(1) > x(end);  if flip_x, x = fliplr(x); end
flip_y = y(1) > y(end);  if flip_y, y = flipud(y); end

% --- potong area sekitar titik + margin (margin > radius kernel) ---
rad    = ceil(4*sigma);
margin = rad + 10;
dxd = mean(diff(x));  dyd = mean(diff(y));
i0 = max(1,        find(x >= min(lon_q) - 1e-9, 1, 'first') - margin);
i1 = min(numel(x), find(x <= max(lon_q) + 1e-9, 1, 'last')  + margin);
j0 = max(1,        find(y >= min(lat_q) - 1e-9, 1, 'first') - margin);
j1 = min(numel(y), find(y <= max(lat_q) + 1e-9, 1, 'last')  + margin);

% indeks asli di file (sebelum dibalik)
if flip_x, ri = numel(x) - [i1 i0] + 1; else, ri = [i0 i1]; end
if flip_y, rj = numel(y) - [j1 j0] + 1; else, rj = [j0 j1]; end
zinfo = info.Variables(strcmp(vn, nz_));
dimNames = {zinfo.Dimensions.Name};
if strcmp(dimNames{1}, nx_) || any(strcmp(dimNames{1}, {'x','lon','longitude'}))
    Z = double(ncread(dem_file, nz_, [ri(1) rj(1)], [diff(ri)+1, diff(rj)+1]))';   % [lat x lon]
else
    Z = double(ncread(dem_file, nz_, [rj(1) ri(1)], [diff(rj)+1, diff(ri)+1]));
end
if flip_x, Z = fliplr(Z); end
if flip_y, Z = flipud(Z); end
xs = x(i0:i1);  ys = y(j0:j1);
if any(isnan(Z(:)))
    warning('hitung_slope_aspect: %d sel DEM bernilai NaN diisi nilai tetangga.', nnz(isnan(Z)));
    Z = isiNaN(Z);
end

% --- smoothing Gaussian (separabel) ---
if sigma > 0
    t = -rad:rad;
    g = exp(-0.5*(t/sigma).^2);  g = g / sum(g);
    Zp = [repmat(Z(:,1),1,rad), Z, repmat(Z(:,end),1,rad)];
    Zp = [repmat(Zp(1,:),rad,1); Zp; repmat(Zp(end,:),rad,1)];
    Zs = conv2(g, g, Zp, 'valid');
else
    Zs = Z;
end

% --- gradien dalam meter ---
dx_m = dxd * 111320 * cosd(lat0);
dy_m = dyd * 111320;
[gE, gN] = gradient(Zs, dx_m, dy_m);        % gE = dz/dE, gN = dz/dN

sE = interp2(xs, ys, gE, lon_q, lat_q, 'linear');
sN = interp2(xs, ys, gN, lon_q, lat_q, 'linear');
slope  = atand(hypot(sE, sN));
aspect = mod(atan2d(-sE, -sN), 360);
end

function nama = cariNama(daftar, calon)
nama = '';
for k = 1:numel(calon)
    if any(strcmp(daftar, calon{k})), nama = calon{k}; return; end
end
end

function Z = isiNaN(Z)
% isi NaN dengan rata-rata tetangga secara berulang
while any(isnan(Z(:)))
    M = ~isnan(Z);  Z0 = Z;  Z0(~M) = 0;
    k = [0 1 0; 1 0 1; 0 1 0];
    s = conv2(Z0, k, 'same');  n = conv2(double(M), k, 'same');
    isi = ~M & n > 0;
    if ~any(isi(:)), break; end
    Z(isi) = s(isi) ./ n(isi);
end
end
