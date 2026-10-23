function [hasil] = SPF_aspect_joughin_vce(value_type, varargin)
% Dekomposisi 3D (dU, dE, dN) dari InSAR ascending + descending dengan kendala
% surface-parallel flow (SPF, Joughin et al. 1998) dan slope aspect
% (Isya et al. 2019), bobot VCE, serta opsi pengisian piksel ill-conditioned
% (regularisasi Tikhonov atau GMT surface).
%
% 12.02.2019 NI ; pseudo 3D dari LOS asc, dsc dan SPF
% 05.07.2019 NI ; kontrol SPF dengan slope aspect
% 15.07.2019 NI ; kontrol optimal dengan VCE
% 24.10.2019 NI ; dengan/tanpa sensitivity (R-index)
% Revisi        ; SPF pers. (18) Joughin dengan gradien lereng, Joughin terkoreksi,
%                 baris aspect benar, VCE Helmert/Förstner global, Tikhonov /
%                 GMT surface, dan mode TIME SERIES (displacement kumulatif).
%
% PEMAKAIAN:
%   SPF_aspect_joughin_vce('mean') : data_match.mat berisi var_vector (npix x 2)
%   Pengaturan bisa ditimpa: SPF_aspect_joughin_vce('ts','metode','vce_tikhonov','P_aspect',0.5)
%   SPF_aspect_joughin_vce('ts')   : data_match.mat berisi var_vector_asc/_dsc
%                                    (npix x nt) dan interpolate.mat berisi range
%                                    (tanggal, datenum). Hasil: displacement
%                                    kumulatif terhadap epoch pertama.
% Butuh: vce.m, joughin_spf.m, dekomposisi_tikhonov.m ('vce_tikhonov'),
%        isi_gmt_surface.m + GMT ('vce_surface'), hitung_slope_aspect.m ('dem').

hasil = [];
if ~(strcmp(value_type,'mean') || strcmp(value_type,'ts'))
    disp('Step is unknown (pakai ''mean'' atau ''ts'')');  return
end

%% ======================= PENGATURAN =======================
sumber_lereng = 'dem';    % 'auto' : var_slope + var_angle(:,5) dari data_match.mat bila ada,
                           %          selain itu dari DEM
                           % 'data' : var_slope/var_angle(:,5) (data_match) atau data.mat
                           % 'dem'  : dihitung ulang dari dem.grd
dem_file      = 'dem.grd';
pakai_aspect  = true;      % false untuk gletser / aliran yang tidak mengikuti aspect lokal
metode        = 'vce_tikhonov';     % 'vce'          : LOS + SPF + aspect, bobot VCE (per piksel)
                           % 'vce_tikhonov' : bobot VCE + regularisasi spasial Laplacian
                           % 'vce_surface'  : hasil 'vce', piksel ill-conditioned diisi GMT surface
                           % 'joughin'      : Joughin (1998) murni (LOS asc+dsc + SPF)
vce_epoch     = 'semua';   % time series: 'akhir' = bobot VCE diestimasi sekali dari epoch
                           %   terakhir (sinyal kumulatif terbesar) lalu dipakai semua epoch;
                           %   'semua' = VCE per epoch (lebih lambat)
ref_epoch1    = true;      % time series: kurangi epoch pertama -> kumulatif mulai dari 0
lambda_tik    = 1.0;       % kekuatan regularisasi Tikhonov. Bila diberi VEKTOR, mis.
                           % logspace(-1,1.5,12), lambda dipilih otomatis dari sudut
                           % kurva-L pada epoch terakhir (gambar kurva_L_tikhonov.png)
tension_gmt   = 0.35;      % faktor tegangan GMT surface
R_min         = 0.05;      % bobot minimum LOS (R-index <= R_min)
P_spf         = 0.41;         % bobot awal grup SPF (pers. 14 naskah: 0.41). Diestimasi ulang oleh VCE
P_aspect      = 0.155;         % bobot grup aspect   (pers. 14 naskah: 0.155). Grup ini DI-FIX di VCE,
                           % sehingga nilai ini tetap dipakai sampai solusi akhir
cond_lim      = 10;        % penanda geometri lemah
maxiter       = 300;
% kontrol kualitas (QC): piksel yang ditandai di-NaN-kan pada dX_ts_new (dX_ts tetap utuh)
qc_aktif      = true;      % false = QC OFF: tidak ada piksel yang di-NaN-kan (penanda QC tetap
                           %         dihitung & disimpan sebagai informasi: flag_C, flag_naik, ...)
mask_qc       = {'C','naik'};  % 'C'   : LOS asc & dsc hampir buta terhadap gerak searah lereng
                             % 'naik': gerak horizontal berlawanan arah turun lereng (khusus longsor)
                             % 'res' : (opsional) residual model besar. Hati-hati: deformasi nyata
                             %         yang tidak memenuhi SPF (mis. amblesan di zona depletion)
                             %         juga memberi residual besar, sehingga ikut termasker
C_min         = 0.2;       % batas sensitivitas |C| (Notti et al. 2014), max(|C_asc|,|C_dsc|) < C_min
z_res         = 3;         % batas residual ternormalisasi (median + z_res*NMAD)
estimator     = 'helmert'; % 'helmert' (otomatis ke Förstner bila varians <= 0) / 'forstner'

% ---- pengaturan dapat ditimpa saat memanggil fungsi (pasangan nama-nilai), mis.
%      SPF_aspect_joughin_vce('ts', 'metode','vce_tikhonov', 'P_aspect',0.5)
if mod(numel(varargin), 2) ~= 0
    error('Argumen tambahan harus berpasangan: nama, nilai.');
end
for a_ = 1:2:numel(varargin)
    nama_ = varargin{a_};
    if ~ischar(nama_) || ~exist(nama_, 'var') || any(strcmp(nama_, {'value_type','varargin','hasil'}))
        error('Pengaturan "%s" tidak dikenal.', num2str(nama_));
    end
    eval([nama_ ' = varargin{a_+1};']);
    if ischar(varargin{a_+1}), fprintf('  %s = %s\n', nama_, varargin{a_+1});
    else, fprintf('  %s = %s\n', nama_, mat2str(varargin{a_+1})); end
end
clear a_ nama_

vce_parm   = strtrim(fileread('vce_parm.txt'));
value_parm = dlmread('value_parm.txt');
value_convergent = value_parm(2);

%% ======================= DATA =======================
D = load('data_match.mat');
var_angle = D.var_angle;  var_lonlat = D.var_lonlat;
if strcmp(value_type, 'ts')
    disp('Dekomposisi 3D TIME SERIES (displacement kumulatif)')
    if isfield(D, 'uw_correct_asc')
        LOSa = D.uw_correct_asc;  LOSd = D.uw_correct_dsc;
    else
        LOSa = D.var_vector_asc;  LOSd = D.var_vector_dsc;
    end
    I = load('interpolate.mat', 'range');
    tanggal = I.range(:);
    if size(LOSa,2) ~= numel(tanggal) || size(LOSd,2) ~= numel(tanggal)
        error('Jumlah kolom var_vector_asc/dsc (%d/%d) tidak sama dengan range (%d).', ...
              size(LOSa,2), size(LOSd,2), numel(tanggal));
    end
    if ref_epoch1
        LOSa = LOSa - LOSa(:,1);  LOSd = LOSd - LOSd(:,1);
    end
else
    disp('Dekomposisi 3D kecepatan / displacement rata-rata')
    LOSa = D.var_vector(:,1);  LOSd = D.var_vector(:,2);
    tanggal = 1;
end
LOSa = double(LOSa);  LOSd = double(LOSd);
[npix, nt] = size(LOSa);
fprintf('%d piksel, %d epoch, metode %s\n', npix, nt, metode);

% slope & aspect -> gradien lereng
if strcmp(sumber_lereng, 'auto')
    if isfield(D, 'var_slope'), sumber_lereng = 'data'; else, sumber_lereng = 'dem'; end
end
if strcmp(sumber_lereng, 'dem')
    [sE, sN, slope, aspect] = hitung_slope_aspect(dem_file, var_lonlat(:,1), var_lonlat(:,2), 2);
else
    if isfield(D, 'var_slope')
        slope = abs(D.var_slope(:));  aspect = var_angle(:,5);
    else
        T = load('data.mat', 'slope', 'aspect');
        slope = T.slope(:,3);  aspect = T.aspect(:,3);
    end
    sE = -tand(slope).*sind(aspect);        % dz/dE
    sN = -tand(slope).*cosd(aspect);        % dz/dN
end
fprintf('Slope/aspect dari: %s\n', sumber_lereng);
aspect = mod(aspect, 360);                       % mis. nilai -8 atau 374 -> 352 / 14
flag_lereng = isnan(slope) | isnan(aspect) | isnan(sE) | isnan(sN);
if any(flag_lereng)
    % piksel datar / void DEM: aspect tidak terdefinisi. Tanpa penanganan, satu NaN
    % membuat seluruh solusi Tikhonov (satu sistem untuk semua piksel) menjadi NaN.
    warning(['%d piksel (%.1f%%) tanpa slope/aspect (permukaan datar atau void DEM): ' ...
             'slope = 0 (SPF -> dU = 0)%s.'], nnz(flag_lereng), 100*mean(flag_lereng), ...
            ternary(pakai_aspect, ' dan piksel ini ditandai QC', ''));
    slope(flag_lereng) = 0;  aspect(flag_lereng) = 0;
    sE(flag_lereng) = 0;     sN(flag_lereng) = 0;
end
theta = 90 - aspect;

% bobot awal (sensitivity = R-index, Notti et al. 2014)
P0 = ones(4, npix);
P0(3,:) = P_spf;  P0(4,:) = P_aspect;
if strcmp(vce_parm, 'sensitivity')
    Ra = [];  Rd = [];
    if isfield(D, 'rindex_asc')
        Ra = D.rindex_asc(:);  Rd = D.rindex_dsc(:);  sumber_R = 'data_match.mat';
    elseif exist('V_slope.mat', 'file')
        R = load('V_slope.mat', 'r_index_asc', 'r_index_dsc');
        Ra = R.r_index_asc(:);  Rd = R.r_index_dsc(:);  sumber_R = 'V_slope.mat';
    end
    if numel(Ra) ~= npix || numel(Rd) ~= npix
        % R-index tidak ada atau milik data lain (jumlah piksel berbeda):
        % hitung ulang dari geometri piksel ini (rumus sama dengan vslope.m, Notti et al. 2014)
        if ~isempty(Ra)
            warning(['R-index di %s berisi %d piksel, data ini %d piksel (kemungkinan dari ' ...
                     'area/studi lain). R-index dihitung ulang dari slope, aspect dan geometri.'], ...
                    sumber_R, numel(Ra), npix);
        end
        Ra = hitung_rindex(slope, aspect, var_angle(:,1), var_angle(:,3));
        Rd = hitung_rindex(slope, aspect, var_angle(:,2), var_angle(:,4));
    end
    P0(1,:) = max(Ra, R_min)';  P0(2,:) = max(Rd, R_min)';
    fprintf('R-index <= %.2f: %d piksel asc, %d piksel dsc\n', R_min, nnz(Ra <= R_min), nnz(Rd <= R_min));
end

%% ======================= MATRIKS DESAIN =======================
% parameter [dU dE dN]; baris: LOS asc, LOS dsc (positif ke satelit),
% SPF (-dU + sE dE + sN dN = 0), aspect (sin(theta) dE - cos(theta) dN = 0)
ia = var_angle(:,3);  id = var_angle(:,4);  ha = var_angle(:,1);  hd = var_angle(:,2);
Ball = zeros(4,3,npix);
Ball(1,:,:) = [cosd(ia), -sind(ia).*sind(ha+90), -sind(ia).*cosd(ha+90)]';
Ball(2,:,:) = [cosd(id), -sind(id).*sind(hd+90), -sind(id).*cosd(hd+90)]';
Ball(3,:,:) = [-ones(npix,1), sE, sN]';
Ball(4,:,:) = [zeros(npix,1), sind(theta), -cosd(theta)]';

if pakai_aspect
    grp = [1 1 2 3];  ref = 1;  fixed = 3;       % aspect di-fix
    baris_cond = [1 2 4];                         % geometri LOS + aspect
else
    Ball = Ball(1:3,:,:);  P0 = P0(1:3,:);
    grp = [1 1 2];    ref = 1;  fixed = 2;       % tanpa redundansi: bobot tidak diestimasi
    baris_cond = [1 2 3];                         % geometri LOS + SPF
end
nobs = size(Ball,1);
kond = zeros(npix,1);
for c = 1:npix, kond(c) = cond(Ball(baris_cond,:,c)); end
flag_ill = kond >= cond_lim;
fprintf('Piksel dengan cond >= %g: %d dari %d\n', cond_lim, nnz(flag_ill), npix);

% observasi per epoch: nobs x npix x nt (baris kendala = 0)
Aall = zeros(nobs, npix, nt);
Aall(1,:,:) = reshape(LOSa, 1, npix, nt);
Aall(2,:,:) = reshape(LOSd, 1, npix, nt);

%% ======================= JOUGHIN (pembanding, selalu dihitung) =======================
[dU_joughin, dE_joughin, dN_joughin] = joughin_spf(LOSa, LOSd, ha, hd, ia, id, sE, sN);

%% ======================= SOLUSI =======================
ep_valid = 1:nt;
if strcmp(value_type,'ts') && ref_epoch1, ep_valid = 2:nt; end   % epoch 1 = 0
dU_ts = zeros(npix,nt);  dE_ts = dU_ts;  dN_ts = dU_ts;
dU_pp = dU_ts;  dE_pp = dU_ts;  dN_pp = dU_ts;                    % per piksel (sebelum pengisian)
k_ts  = nan(max(grp), nt);
vce_sta = {};

if strcmp(metode, 'joughin')
    disp('Metode joughin: solusi langsung (Joughin et al. 1998), VCE tidak dijalankan -> tanpa iterasi.')
    dU_ts = dU_joughin;  dE_ts = dE_joughin;  dN_ts = dN_joughin;
    dU_pp = dU_ts;  dE_pp = dE_ts;  dN_pp = dN_ts;
else
    % --- bobot VCE ---
    if strcmp(vce_epoch, 'semua'), ep_vce = ep_valid; else, ep_vce = ep_valid(end); end
    for n = ep_vce
        if isempty(setdiff(1:max(grp), [ref fixed]))
            disp(['VCE tidak dijalankan: tanpa baris aspect sistem LOS + SPF tidak memiliki ' ...
                  'redundansi (3 persamaan, 3 parameter), bobot tidak dapat diestimasi -> k = 1.'])
            k = ones(max(grp),1);  st = [];
        else
            fprintf('--- VCE epoch %d/%d ---\n', n, nt);
            [~, k, st] = vce(Ball, Aall(:,:,n), P0, grp, ref, fixed, value_convergent, maxiter, estimator);
        end
        k_ts(:,n) = k;  vce_sta{end+1} = st; %#ok<AGROW>
    end
    if ~strcmp(vce_epoch, 'semua'), k_ts(:, ep_valid) = repmat(k, 1, numel(ep_valid)); end
    fprintf('Faktor bobot akhir: %s\n', num2str(k_ts(:, ep_valid(end))', '%10.4f'));

    % --- solusi per piksel ---
    for n = ep_valid
        x = ls_perpiksel(Ball, Aall(:,:,n), k_ts(grp,n) .* P0);
        dU_pp(:,n) = x(1,:)';  dE_pp(:,n) = x(2,:)';  dN_pp(:,n) = x(3,:)';
    end
    dU_ts = dU_pp;  dE_ts = dE_pp;  dN_ts = dN_pp;

    % --- pengisian piksel ill-conditioned ---
    if strcmp(metode, 'vce_tikhonov')
        if numel(lambda_tik) > 1
            lambda_tik = pilih_lambda(Ball, Aall(:,:,ep_valid(end)), k_ts(grp,ep_valid(end)).*P0, ...
                                      var_lonlat, lambda_tik);
        end
        if strcmp(vce_epoch, 'semua')
            for n = ep_valid
                x = dekomposisi_tikhonov(Ball, Aall(:,:,n), k_ts(grp,n).*P0, ...
                                         var_lonlat(:,1), var_lonlat(:,2), lambda_tik);
                dU_ts(:,n) = x(1,:)';  dE_ts(:,n) = x(2,:)';  dN_ts(:,n) = x(3,:)';
            end
        else
            % bobot sama untuk semua epoch -> satu faktorisasi, semua epoch sekaligus
            [x, info_tik] = dekomposisi_tikhonov(Ball, Aall(:,:,ep_valid), k_ts(grp,ep_valid(end)).*P0, ...
                                                 var_lonlat(:,1), var_lonlat(:,2), lambda_tik);
            dU_ts(:,ep_valid) = reshape(x(1,:,:), npix, []);
            dE_ts(:,ep_valid) = reshape(x(2,:,:), npix, []);
            dN_ts(:,ep_valid) = reshape(x(3,:,:), npix, []);
            fprintf('Tikhonov lambda = %g (epoch terakhir: residu %.3f, kekasaran %.3f)\n', ...
                    lambda_tik, info_tik.residu(end), info_tik.kasar(end));
        end
    elseif strcmp(metode, 'vce_surface')
        ne = numel(ep_valid);
        Xf = isi_gmt_surface([dU_pp(:,ep_valid) dE_pp(:,ep_valid) dN_pp(:,ep_valid)], ...
                             var_lonlat(:,1), var_lonlat(:,2), flag_ill, tension_gmt);
        dU_ts(:,ep_valid) = Xf(:, 1:ne);
        dE_ts(:,ep_valid) = Xf(:, ne+1:2*ne);
        dN_ts(:,ep_valid) = Xf(:, 2*ne+1:3*ne);
        fprintf('%d piksel ill-conditioned diisi GMT surface (T = %g) di %d epoch\n', ...
                nnz(flag_ill), tension_gmt, ne);
    end
end

%% ======================= KONTROL KUALITAS =======================
% (1) sensitivitas LOS terhadap gerak searah lereng: C = L . u_turun (Notti et al. 2014)
u_t = [cosd(slope).*sind(aspect), cosd(slope).*cosd(aspect), -sind(slope)];   % [E N U]
L_a = [-sind(ia).*sind(ha+90), -sind(ia).*cosd(ha+90), cosd(ia)];
L_d = [-sind(id).*sind(hd+90), -sind(id).*cosd(hd+90), cosd(id)];
C_asc = sum(L_a.*u_t, 2);  C_dsc = sum(L_d.*u_t, 2);
flag_C = max(abs(C_asc), abs(C_dsc)) < C_min;
% (2) konsistensi data: residual berbobot solusi per piksel di epoch terakhir
ne_ = ep_valid(end);
if strcmp(metode, 'joughin') || all(isnan(k_ts(:,ne_))), kw = ones(max(grp),1); else, kw = k_ts(:,ne_); end
Wq = kw(grp(:)) .* P0;
xq = [dU_pp(:,ne_) dE_pp(:,ne_) dN_pp(:,ne_)]';
vq = squeeze(sum(Ball .* reshape(xq, 1, 3, npix), 2)) - Aall(:,:,ne_);
res_qc = sqrt(sum(Wq .* vq.^2, 1))';
nmad = 1.4826*median(abs(res_qc - median(res_qc)));
z_qc = (res_qc - median(res_qc)) / max(nmad, eps);
flag_res = z_qc > z_res;
% (3) gerak naik lereng (tidak wajar untuk longsor)
Hq = hypot(dE_ts(:,ne_), dN_ts(:,ne_));
arah = (dE_ts(:,ne_).*sind(aspect) + dN_ts(:,ne_).*cosd(aspect)) ./ max(Hq, eps);
flag_naik = arah < -0.5 & Hq > 3*1.4826*median(abs(Hq - median(Hq)));
fprintf('QC: |C| < %.2f: %d piksel | residual z > %g: %d piksel | gerak naik lereng: %d piksel\n', ...
        C_min, nnz(flag_C), z_res, nnz(flag_res), nnz(flag_naik));
flag_qc = false(npix,1);
if ~qc_aktif
    mask_qc = {};                                    % QC dimatikan
end
if pakai_aspect, flag_qc = flag_qc | flag_lereng; end   % arah aspect tidak terdefinisi
if ~pakai_aspect && any(strcmp(mask_qc, 'naik'))
    mask_qc = setdiff(mask_qc, {'naik'});            % gletser: arah aliran tidak harus = aspect lokal
    disp('QC: kriteria ''naik'' dinonaktifkan karena pakai_aspect = false');
end
if any(strcmp(mask_qc, 'C')),    flag_qc = flag_qc | flag_C;    end
if any(strcmp(mask_qc, 'res')),  flag_qc = flag_qc | flag_res;  end
if any(strcmp(mask_qc, 'naik')), flag_qc = flag_qc | flag_naik; end
if qc_aktif
    fprintf('QC: %d piksel (%.1f%%) di-NaN-kan pada dX_ts_new\n', nnz(flag_qc), 100*mean(flag_qc));
else
    disp('QC OFF: semua piksel dipertahankan (dX_ts_new = dX_ts); penanda QC hanya disimpan.')
end

%% ======================= SIMPAN =======================
dU_ts_new = dU_ts;  dE_ts_new = dE_ts;  dN_ts_new = dN_ts;     % tanpa koreksi kuadran
dU_ts_new(flag_qc,:) = NaN;  dE_ts_new(flag_qc,:) = NaN;  dN_ts_new(flag_qc,:) = NaN;
m_vce_piksel = [dU_pp(:,end)'; dE_pp(:,end)'; dN_pp(:,end)'];  % epoch terakhir (untuk plot_compare)
if strcmp(value_type, 'ts')
    fid = fopen('date.in', 'w');
    ds = cellstr(datestr(tanggal));
    fprintf(fid, '%s\n', ds{:});
    fclose(fid);
end
vv = dU_ts_new(~isnan(dU_ts_new));  scale = [min(vv) max(vv)];
dlmwrite('dU_ts_3d_SF_opt.txt', dU_ts_new, 'precision', 8, 'delimiter', ' ');
dlmwrite('dE_ts_3d_SF_opt.txt', dE_ts_new, 'precision', 8, 'delimiter', ' ');
dlmwrite('dN_ts_3d_SF_opt.txt', dN_ts_new, 'precision', 8, 'delimiter', ' ');
dlmwrite('lonlat.txt', var_lonlat, 'precision', 10, 'delimiter', ' ');
dlmwrite('ver_scale_3d_SF_opt.txt', scale, 'precision', 8, 'delimiter', ' ');
save('generate_3d_SF_opt.mat', 'dU_ts', 'dU_ts_new', 'dE_ts', 'dE_ts_new', 'dN_ts', 'dN_ts_new', ...
     'dU_pp', 'dE_pp', 'dN_pp', 'dU_joughin', 'dE_joughin', 'dN_joughin', 'kond', 'flag_ill', ...
     'slope', 'aspect', 'sE', 'sN', 'metode', 'm_vce_piksel', 'tanggal', 'k_ts', 'var_lonlat', ...
     'flag_qc', 'flag_C', 'flag_res', 'flag_naik', 'flag_lereng', 'C_asc', 'C_dsc', 'z_qc');
save('vce_SF_opt.mat', 'vce_sta', 'k_ts', 'grp');
if nargout > 0, hasil = cat(3, dU_ts, dE_ts, dN_ts); end
disp('Selesai: generate_3d_SF_opt.mat dan vce_SF_opt.mat tersimpan.')
end

%% ======================= FUNGSI LOKAL =======================
function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function R = hitung_rindex(slope, aspect, heading, inc)
% R-index (Notti et al. 2014), identik dengan vslope.m
los_azi = mod(heading(:) + 90, 360);
R = -sind(slope(:) .* sind(aspect(:) + 270 - los_azi) - inc(:));
end

function lam = pilih_lambda(Ball, A, W, ll, daftar)
% kurva-L: norma residual vs norma kekasaran untuk beberapa lambda,
% lambda terpilih = titik dengan kelengkungan maksimum (sudut kurva)
daftar = sort(daftar(:))';
rs = zeros(size(daftar));  ks = rs;
for j = 1:numel(daftar)
    [~, info] = dekomposisi_tikhonov(Ball, A, W, ll(:,1), ll(:,2), daftar(j));
    rs(j) = info.residu;  ks(j) = info.kasar;
    fprintf('  kurva-L: lambda = %8.3f | residu %10.3f | kekasaran %10.3f\n', daftar(j), rs(j), ks(j));
end
x = log10(rs);  y = log10(ks);
kel = zeros(size(x));
for j = 2:numel(x)-1                      % kelengkungan Menger dari tiga titik berurutan
    a = [x(j-1) y(j-1)];  b = [x(j) y(j)];  c = [x(j+1) y(j+1)];
    luas2 = (b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1));
    kel(j) = 2*luas2 / (norm(b-a)*norm(c-b)*norm(c-a) + eps);
end
[~, jm] = max(abs(kel));
lam = daftar(jm);
fprintf('  lambda terpilih (sudut kurva-L): %.3f\n', lam);
f = figure('Color','w');
loglog(rs, ks, 'o-'); hold on; loglog(rs(jm), ks(jm), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
for j = 1:numel(daftar), text(rs(j), ks(j), sprintf('  %.2g', daftar(j)), 'FontSize', 8); end
grid on; xlabel('Residual norm ||W^{1/2}(Ax - y)||'); ylabel('Roughness norm ||Lx||');
title(sprintf('L-curve (selected \\lambda = %.3g)', lam));
print(f, '-dpng', '-r200', 'kurva_L_tikhonov.png');
end

function X = ls_perpiksel(Ball, A, W)
% kuadrat terkecil berbobot per piksel (vektor), A dan W: nobs x npix
[nobs, ~, npix] = size(Ball);
b = cell(1,3);
for a = 1:3, b{a} = reshape(Ball(:,a,:), nobs, npix); end
N11 = sum(W.*b{1}.*b{1},1); N12 = sum(W.*b{1}.*b{2},1); N13 = sum(W.*b{1}.*b{3},1);
N22 = sum(W.*b{2}.*b{2},1); N23 = sum(W.*b{2}.*b{3},1); N33 = sum(W.*b{3}.*b{3},1);
r1 = sum(W.*b{1}.*A,1);  r2 = sum(W.*b{2}.*A,1);  r3 = sum(W.*b{3}.*A,1);
C11 = N22.*N33 - N23.^2;  C12 = N13.*N23 - N12.*N33;  C13 = N12.*N23 - N13.*N22;
C22 = N11.*N33 - N13.^2;  C23 = N12.*N13 - N11.*N23;  C33 = N11.*N22 - N12.^2;
Dt  = N11.*C11 + N12.*C12 + N13.*C13;
X = [ (C11.*r1 + C12.*r2 + C13.*r3) ./ Dt;
      (C12.*r1 + C22.*r2 + C23.*r3) ./ Dt;
      (C13.*r1 + C23.*r2 + C33.*r3) ./ Dt ];
end
