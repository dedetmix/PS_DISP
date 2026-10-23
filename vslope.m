function [asc, dsc] = vslope(value_type, varargin)
% VSLOPE  R-index, faktor C dan kecepatan sepanjang lereng (V_slope)
%
% 11.06.2018 NI ; R-index dan Vslope berdasarkan C
% Revisi        ; fungsi bersama asc/dsc, penanda layover/shadow (R-index <= 0),
%                 penyaringan |C| kecil agar V_slope tidak meledak,
%                 konvensi tanda didokumentasikan, pemuatan data yang tidak
%                 dipakai dihapus.
%
% Sumber: Notti et al. (2014), "A methodology for improving landslide PSI data
%         analysis", Int. J. Remote Sensing.
%         Pizarro et al. (2017), MDPI Remote Sensing.
%
% PEMAKAIAN:
%   vslope('asc'), vslope('dsc')                 : R-index, C, V_slope (kecepatan rata-rata)
%   vslope('asc_ts_surf'), vslope('dsc_ts_surf') : V_slope time series (Surface)
%   vslope('asc_ts_nn'),   vslope('dsc_ts_nn')   : V_slope time series (Near-Neighbour)
%   Argumen kedua opsional: batas |C| minimum, mis. vslope('asc', 0.3). Default 0.2.
%
% KONVENSI (sudah dicek dengan geometri eksak):
%   LOS positif = mendekati satelit.
%   C = -(L . u), L = vektor satuan LOS ke satelit, u = vektor satuan TURUN lereng.
%   V_slope = LOS / C  ->  NEGATIF untuk gerak menuruni lereng (konvensi Notti).
%   R-index ~ sin(sudut datang lokal); R <= 0 berarti layover / tidak terlihat.
%
% OUTPUT (ditambahkan ke V_slope.mat):
%   r_index_xxx, C_xxx, vslope_xxx     : seperti versi lama
%   layover_xxx                        : true bila R-index <= 0
%   validC_xxx                         : true bila |C| >= C_min (V_slope = NaN bila false)
%   xxx_ts_surf / xxx_ts_nn            : V_slope time series (NaN bila |C| < C_min)

asc = [];  dsc = [];
C_min = 0.2;
if nargin > 1 && ~isempty(varargin{1}), C_min = varargin{1}; end

switch value_type
    case {'asc', 'dsc'}
        trk = value_type;
        fprintf('V slope dari %s\n', trk);
        D = load('data.mat', 'aspect', 'slope', ['az_angle_' trk], ['inc_angle_' trk], ['los_' trk]);
        azi = D.(['az_angle_' trk])(:,3);
        inc = D.(['inc_angle_' trk])(:,3);
        los = D.(['los_' trk])(:,3);
        [R, C] = hitung_R_C(D.slope(:,3), D.aspect(:,3), azi, inc);

        validC = abs(C) >= C_min;
        vs = los ./ C;
        vs(~validC) = NaN;
        layover = R <= 0;

        fprintf('  R-index <= 0 (layover/shadow): %d piksel\n', nnz(layover));
        fprintf('  |C| < %.2f (V_slope tidak dihitung): %d piksel (%.1f%%)\n', ...
                C_min, nnz(~validC), 100*mean(~validC));

        S.(['r_index_' trk]) = R;
        S.(['C_' trk])       = C;
        S.(['vslope_' trk])  = vs;
        S.(['layover_' trk]) = layover;
        S.(['validC_' trk])  = validC;
        simpan('V_slope.mat', S);
        if strcmp(trk, 'asc'), asc = vs; else, dsc = vs; end

    case {'asc_ts_surf', 'dsc_ts_surf', 'asc_ts_nn', 'dsc_ts_nn'}
        % jalankan di folder time series; vslope('asc'/'dsc') harus dihitung dulu
        trk  = value_type(1:3);
        jenis = value_type(8:end);                     % 'surf' atau 'nn'
        fprintf('V slope time series (%s) dari %s\n', jenis, trk);
        I = load('interpolate.mat', [trk '_' jenis], 'range');
        V = load(fullfile('..', 'V_slope.mat'), ['C_' trk]);
        C = V.(['C_' trk]);
        validC = abs(C) >= C_min;
        ts = I.([trk '_' jenis]) ./ C;                  % bagi setiap kolom (tanggal) dengan C
        ts(~validC, :) = NaN;
        fprintf('  %d tanggal, %d piksel dengan |C| < %.2f diberi NaN\n', ...
                numel(I.range), nnz(~validC), C_min);

        S.([trk '_ts_' jenis]) = ts;
        simpan('V_slope.mat', S);
        if strcmp(trk, 'asc'), asc = ts; else, dsc = ts; end

    otherwise
        disp('Please specify asc, dsc, asc_ts_surf, dsc_ts_surf, asc_ts_nn or dsc_ts_nn');
end
end

%% ======================= FUNGSI LOKAL =======================
function [R, C] = hitung_R_C(slope, aspect, azi, inc)
% R-index (Notti et al., 2014) dan faktor C untuk satu orbit.
%   azi : heading satelit (derajat dari utara), inc : sudut datang (derajat)
los_azi = mod(azi + 90, 360);          % azimut arah pandang (right-looking)

% R-index
So = sind(aspect + 270 - los_azi);
R  = -sind(slope .* So - inc);

% komponen vektor satuan LOS (ke arah satelit)
H = cosd(inc);
N = cosd(90 - inc) .* cosd(180 - los_azi);
E = cosd(90 - inc) .* cosd(270 - los_azi);

% C = proyeksi LOS pada arah NAIK lereng = -(L . u_turun)
C = cosd(slope) .* sind(aspect - 90) .* N ...
  - cosd(slope) .* cosd(aspect - 90) .* E ...
  + sind(slope) .* H;
end

function simpan(file, S)
% simpan semua field struct S ke file .mat (tambahkan bila file sudah ada)
if exist(file, 'file')
    save(file, '-struct', 'S', '-append');
else
    save(file, '-struct', 'S');
end
end
