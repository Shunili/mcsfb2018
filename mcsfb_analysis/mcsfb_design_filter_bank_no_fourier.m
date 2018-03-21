function [ filter_bank, shifted_ends] = mcsfb_design_filter_bank_no_fourier( G, num_bands, band_ends, param)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(param,'band_structure')
    param.band_structure = 0; %1: minima of cdf; 0: search in the interval
end

if ~isfield(G,'spectrum_cdf_approx')
    G=gsp_spectrum_cdf_approx(G);
end

if ~isfield(G,'spectrum_pdf_approx')
% compute pdf
    xx = 0:0.001:G.lmax;
    delta=.001;
    G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative
end

filter_bank=cell(num_bands,1);
shifted_ends = zeros(num_bands+1,1);
shifted_ends(1) = 0;
shifted_ends(num_bands+1) = G.lmax;

if param.band_structure %minima of cdf
    % compute pdf minima
    xx = 0:0.001:G.lmax;
    yy = G.spectrum_cdf_approx(xx);
    Dfun = diff(yy)/0.001;% first derivative
    inverted = - Dfun;
    [~, idx] = findpeaks(inverted); % find another way to do it
    pdf_minima = idx*0.001;

    % find the closest values for band_ends
    for i = 1:(length(band_ends)-2)
        [shifted_ends(1+i), idx] = binary_search(band_ends(i+1), pdf_minima);
        %pdf_minima = pdf_minima((idx+1):length(pdf_minima));
    end

    for i = 1:num_bands
        filter_bank{i}=@(x) ((shifted_ends(i) <= x) & (x <= shifted_ends(i+1)));
    end

else %search in the interval
    % search the minima in range 
    for k = 2:num_bands
        shifted_ends(k) = fminbnd(G.spectrum_pdf_approx,(band_ends(k)+band_ends(k-1))/2, (band_ends(k+1)+band_ends(k))/2);
    end

    for l = 1:num_bands
        filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
    end
end
end

