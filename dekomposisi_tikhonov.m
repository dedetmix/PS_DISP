function [m, info] = dekomposisi_tikhonov(Ball, Aall, W, lon, lat, lambda)
% DEKOMPOSISI_TIKHONOV  Dekomposisi 3D seluruh piksel sekaligus dengan
% regularisasi spasial Tikhonov (Laplacian), untuk mengatasi ill-posedness
% di piksel yang geometrinya lemah (arah gerak ~ utara-selatan).
%
% Meminimumkan, untuk setiap epoch:
%   sum_i w_i (b_i*x_p - a_i)^2  +  lambda^2 * ( ||L dU||^2 + ||L dE||^2 + ||L dN||^2 )
% L = Laplacian graf pada piksel yang ADA (tetangga 4-arah). Piksel kosong /
% di luar data tidak dianggap bernilai nol (batas Neumann), sehingga cocok
% untuk data PS/grid yang tidak lengkap.
%
% Input:
%   Ball   : nobs x 3 x npix       matriks desain per piksel, parameter [dU dE dN]
%   Aall   : nobs x npix [x nt]    observasi (nt epoch time series diselesaikan sekaligus)
%   W      : nobs x npix           bobot observasi (mis. k(grp).*P0 dari VCE)
%   lon, lat : koordinat piksel (npix x 1) pada grid teratur
%   lambda : kekuatan regularisasi (pilih dengan kurva-L, lihat info)
% Output:
%   m    : 3 x npix [x nt]  [dU; dE; dN]
%   info : .residu, .kasar (per epoch) untuk kurva-L; .n_tetangga per piksel

[nobs, ~, n] = size(Ball);
nt = size(Aall, 3);
lon = lon(:);  lat = lat(:);

% --- indeks grid & tetangga 4-arah di antara piksel yang ada ---
ux = unique(round(lon*1e7));  uy = unique(round(lat*1e7));
[~, ix] = ismember(round(lon*1e7), ux);
[~, iy] = ismember(round(lat*1e7), uy);
peta = sparse(iy, ix, (1:n)', numel(uy), numel(ux));
I = [];  J = [];
for geser = [0 1; 1 0]'
    ok = (iy + geser(1) <= numel(uy)) & (ix + geser(2) <= numel(ux));
    p  = find(ok);
    q  = full(peta(sub2ind(size(peta), iy(ok) + geser(1), ix(ok) + geser(2))));
    ada = q > 0;
    I = [I; p(ada)];  J = [J; q(ada)]; %#ok<AGROW>
end
Adj = sparse([I; J], [J; I], 1, n, n);
deg = full(sum(Adj, 2));
Lp  = Adj - spdiags(deg, 0, n, n);            % Laplacian graf: sum(tetangga) - deg*x_p
L   = blkdiag(Lp, Lp, Lp);

% --- bagian data (berbobot) ---
sw = sqrt(W);
rows = zeros(nobs*3*n,1);  cols = rows;  vals = rows;  c = 0;
for r = 1:nobs
    for a = 1:3
        idx = c + (1:n);
        rows(idx) = (r-1)*n + (1:n);
        cols(idx) = (a-1)*n + (1:n);
        vals(idx) = reshape(Ball(r,a,:), n, 1) .* sw(r,:)';
        c = c + n;
    end
end
A = sparse(rows, cols, vals, nobs*n, 3*n);
Y = reshape(permute(Aall .* sw, [2 1 3]), nobs*n, nt);   % [obs1 semua piksel; obs2 ...] x epoch

% --- solusi persamaan normal (satu faktorisasi untuk semua epoch) ---
Nm = A'*A + lambda^2 * (L'*L);
X  = Nm \ (A'*Y);                                       % 3n x nt
m  = permute(reshape(X, n, 3, nt), [2 1 3]);            % 3 x n x nt

info.residu = sqrt(sum((A*X - Y).^2, 1));
info.kasar  = sqrt(sum((L*X).^2, 1));
info.lambda = lambda;
info.n_tetangga = deg;
end
