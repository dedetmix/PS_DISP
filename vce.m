function [m, k, stat] = vce(Ball, Aall, P0, grp, ref, fixed, tol, maxiter, estimator, redaman, r_min)
% VCE  Variance Component Estimation metode HELMERT (versi terkoreksi dari
% vce.m lama), dihitung GLOBAL untuk semua piksel sekaligus.
%
% Perbaikan terhadap vce.m lama:
%  1. residual tiap grup: v_i = b_i*m - A_mat(i)   (lama: dikurangi seluruh A_mat)
%  2. diagonal S: n_i - 2 tr(N^-1 N_i) + tr(N^-1 N_i N^-1 N_i)
%     (lama: n_i - 2 tr(N^-1 N_i) + tr(N^-1 N_i), lalu dikuadratkan)
%  3. n_i = jumlah observasi grup (lama: 1, 0.41, 0.155 -> bukan jumlah observasi)
%  4. pembaruan bobot relatif terhadap grup acuan: k_i <- k_i*theta_ref/theta_i
%     (lama: pada cabang "else" bobot asc & dsc tertukar)
%  5. kriteria konvergensi relatif: max|theta_i/theta_ref - 1| < tol
%     (lama: max(theta) - min(theta), bergantung satuan)
%  6. dihitung global: per piksel hanya ada 4 observasi untuk 3 parameter
%     (redundansi 1), sehingga 4 komponen varians TIDAK dapat diestimasi per piksel.
%
% Input:
%   Ball    : nobs x 3 x npix  matriks desain tiap piksel ([dU dE dN])
%   Aall    : nobs x npix      observasi
%   P0      : nobs x npix      bobot awal a priori
%   grp     : 1 x nobs         nomor grup tiap baris, mis. [1 1 2 3]
%   ref     : grup acuan;  fixed : grup yang bobotnya tidak diestimasi
%   tol, maxiter : toleransi relatif & iterasi maksimum
%   estimator : 'helmert' (default) atau 'forstner'. Helmert otomatis diganti
%               Förstner pada iterasi di mana hasilnya tidak stabil: varians <= 0,
%               atau berbeda > 5x dari estimasi Förstner (matriks S hampir singular
%               bila grup observasi sulit dipisahkan, mis. LOS vs SPF)
%   redaman   : 0 < a <= 1, pembaruan bobot k <- k .* (theta_ref./theta).^a
%               (default 0.5). a < 1 meredam osilasi; a = 1 = pembaruan penuh (versi lama)
%   r_min     : bila redundansi grup turun di bawah r_min x jumlah observasinya
%               (default 0.02), grup itu praktis menjadi kendala keras: bobotnya
%               tidak dapat diestimasi lagi, sehingga dibekukan pada nilai saat itu
% Output:
%   m    : 3 x npix [dU; dE; dN];  k : faktor bobot per grup
%   stat : [iter err k' theta' r' rmse_LOS MSE std_dU std_dE std_dN] per iterasi
%          (MSE = vTPv/r; std = median antar piksel dari sqrt(MSE*diag(N^-1)))

if nargin < 9 || isempty(estimator), estimator = 'helmert'; end
if nargin < 10 || isempty(redaman), redaman = 0.5; end
if nargin < 11 || isempty(r_min), r_min = 0.02; end
beku = [];
n_ganti = 0;
[nobs, ~, npix] = size(Ball);
grp  = grp(:);
ng   = max(grp);
k    = ones(ng,1);
free = setdiff(1:ng, [ref fixed(:)']);
b    = cell(1,3);
for a = 1:3, b{a} = reshape(Ball(:,a,:), nobs, npix); end
stat = zeros(maxiter, 7 + 3*ng);

for iter = 1:maxiter
    W = k(grp) .* P0;
    % --- matriks normal total N dan per grup N_g (3x3 per piksel) ---
    Ng = zeros(3,3,npix,ng);
    for g = 1:ng
        idx = (grp == g);
        for a = 1:3
            for c = a:3
                Ng(a,c,:,g) = reshape(sum(W(idx,:).*b{a}(idx,:).*b{c}(idx,:), 1), 1,1,npix);
                Ng(c,a,:,g) = Ng(a,c,:,g);
            end
        end
    end
    N = sum(Ng, 4);
    Ni = inv3(N);

    % --- solusi & residual ---
    rhs = zeros(3,npix);
    for a = 1:3, rhs(a,:) = sum(W.*b{a}.*Aall, 1); end
    x = zeros(3,npix);
    for a = 1:3
        x(a,:) = reshape(Ni(a,1,:),1,npix).*rhs(1,:) + reshape(Ni(a,2,:),1,npix).*rhs(2,:) ...
               + reshape(Ni(a,3,:),1,npix).*rhs(3,:);
    end
    m = x;
    v = b{1}.*x(1,:) + b{2}.*x(2,:) + b{3}.*x(3,:) - Aall;

    % --- M_g = N^-1 N_g, lalu matriks S dan vektor w (Helmert) ---
    M = zeros(3,3,npix,ng);
    for g = 1:ng, M(:,:,:,g) = mul3(Ni, Ng(:,:,:,g)); end
    S = zeros(ng);  w = zeros(ng,1);  r = zeros(ng,1);
    for g = 1:ng
        idx  = (grp == g);
        n_g  = nnz(idx) * npix;
        trMg = sum(M(1,1,:,g) + M(2,2,:,g) + M(3,3,:,g));
        r(g) = n_g - trMg;                                   % redundansi grup
        w(g) = sum(sum(W(idx,:).*v(idx,:).^2));              % v' P v
        for h = g:ng
            trMM = sum(reshape(sum(sum(M(:,:,:,g) .* permute(M(:,:,:,h),[2 1 3]),1),2),1,[]));
            if h == g
                S(g,g) = n_g - 2*trMg + trMM;
            else
                S(g,h) = trMM;  S(h,g) = trMM;
            end
        end
    end
    theta_f = w ./ r;                              % Förstner (selalu positif, stabil)
    theta   = theta_f;
    if strcmp(estimator, 'helmert')
        theta_h = S \ w;                          % Helmert: S*theta = w
        rasio   = theta_h ./ theta_f;
        idx_cek = [ref free];
        if all(theta_h(idx_cek) > 0) && all(rasio(idx_cek) > 0.2 & rasio(idx_cek) < 5)
            theta = theta_h;
        else
            n_ganti = n_ganti + 1;                 % Helmert tidak stabil -> pakai Förstner
        end
    end

    % grup dengan redundansi hampir nol -> bekukan (bobot tidak terestimasi)
    for g = free
        if r(g) < r_min * nnz(grp == g) * npix
            beku = [beku g]; %#ok<AGROW>
            fprintf('vce: redundansi grup %d tinggal %.1f%% -> bobot dibekukan (k = %.4g)\n', ...
                    g, 100*r(g)/(nnz(grp == g)*npix), k(g));
        end
    end
    free = setdiff(free, beku);

    if isempty(free), err = 0; else, err = max(abs(theta(free)/theta(ref) - 1)); end
    los = (grp == ref);                            % RMSE residual observasi acuan (LOS)
    rmse_los = sqrt(mean(reshape(v(los,:), [], 1).^2));
    % MSE (variance factor a posteriori, seperti lscov) dan std parameter (median antar piksel)
    mse_tot = sum(W(:) .* v(:).^2) / max(sum(r), eps);
    std_par = zeros(1,3);
    for a = 1:3, std_par(a) = median(sqrt(mse_tot * reshape(Ni(a,a,:), 1, npix))); end
    stat(iter,:) = [iter err k' theta' r' rmse_los mse_tot std_par];
    fprintf('Iterasi %3d | err = %.5f | k = [%s]\n', iter, err, num2str(k','%10.4f'));
    if err < tol, break; end
    upd = (theta(ref) ./ theta) .^ redaman;       % pembaruan teredam
    upd([ref fixed(:)' beku]) = 1;
    k = k .* upd;
end
stat = stat(1:iter,:);
if n_ganti > 0
    fprintf('vce: Helmert tidak stabil pada %d dari %d iterasi, diganti estimasi Förstner.\n', n_ganti, iter);
end
if err >= tol
    warning('VCE belum konvergen setelah %d iterasi (err = %.4g).', maxiter, err);
end
end

function Ni = inv3(N)
% invers 3x3 per piksel (kofaktor), N simetris 3x3xnpix
n11=N(1,1,:); n12=N(1,2,:); n13=N(1,3,:); n22=N(2,2,:); n23=N(2,3,:); n33=N(3,3,:);
c11=n22.*n33-n23.^2; c12=n13.*n23-n12.*n33; c13=n12.*n23-n13.*n22;
c22=n11.*n33-n13.^2; c23=n12.*n13-n11.*n23; c33=n11.*n22-n12.^2;
D = n11.*c11 + n12.*c12 + n13.*c13;
Ni = [c11 c12 c13; c12 c22 c23; c13 c23 c33] ./ D;
end

function C = mul3(A, B)
% perkalian matriks 3x3 per piksel
C = zeros(size(A));
for i = 1:3
    for j = 1:3
        C(i,j,:) = A(i,1,:).*B(1,j,:) + A(i,2,:).*B(2,j,:) + A(i,3,:).*B(3,j,:);
    end
end
end
