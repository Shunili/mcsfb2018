function [filter_bank, shifted_ends, band_ends] = mcsfb_design_filter_bank(G, num_bands,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% return the pdf!
% double check bands!

if ~isfield(param, 'spacing')
    param.spacing = 1; % 1-logarithm; 0-evenly spacing
end

if ~isfield(param, 'spectrum_adapted')
    param.spectrum_adapted = 1;
end

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

filter_bank=cell(num_bands,1);
shifted_ends = zeros(num_bands+1,1);
shifted_ends(1) = 0;
shifted_ends(num_bands+1) = G.lmax+1;

band_ends = zeros(num_bands+1,1);
band_ends(1) = 0;
band_ends(num_bands+1) =G.lmax;

if gsp_check_fourier(G)  
    if param.spacing 
        if param.spectrum_adapted
            % 1.spectrum_adapted=1; spacing = log; check_fourier=1
            for k = 1:num_bands-1
                idx = floor(G.N*(1/2)^k);
                shifted_ends(num_bands-k+1) = (G.e(idx)+G.e(idx+1))/2;
            end
            
            for l = 1:num_bands
                filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
            end
        else
            % 2.spectrum_adapted = 0; spacing = even; check_fourier=1;
            for k = 1:num_bands-1
                eigen = G.lmax*(1/2)^k;
                [~,idx] = min(abs(G.e-eigen));
                shifted_ends(num_bands-k+1) = (G.e(idx)+G.e(idx+1))/2;
            end

            for l = 1:num_bands
                filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x < shifted_ends(l+1)));
            end
        end  
    else
        if param.spectrum_adapted
            % 3. spectrum_adapted=1; spacing = even; check_fourier=1
            for k = 1:num_bands-1   
                idx = floor(G.N/num_bands*k);
                shifted_ends(k+1) = (G.e(idx)+G.e(idx+1))/2;
            end

            for l = 1:num_bands
                filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
            end
        else
            % 4. spectrum_adapted = 0; spacing = even; check_fourier=1; 
            for k = 1:num_bands-1
                eigen = G.lmax/num_bands*k;
                [~,idx] = min(abs(G.e-eigen));
                shifted_ends(k+1) = (G.e(idx)+G.e(idx+1))/2;
            end

            for l = 1:num_bands
                filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
            end
        end
    end
    band_ends = shifted_ends; % for debugging only
else   
%     band_ends = zeros(num_bands+1,1);
%     band_ends(1) = 0;
%     band_ends(num_bands+1) =G.lmax;
    
    if ~isfield(G,'spectrum_cdf_approx')
        [G, ~] = spectral_cdf_approx(G, param);
    end
    
    if ~isfield(G,'spectrum_pdf_approx')
        xx = 0:0.001:G.lmax;
        delta=.1;
        G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative
    end
            
    if param.spacing 
        if param.spectrum_adapted
            % 5.spectrum_adapted=1; spacing = log; check_fourier=0
            % approximate eigenvalues
            for j = 1:num_bands-1
                fun = @(x) (G.spectrum_cdf_approx(x)-1/(2^j));
                band_ends(length(band_ends)-j) = bisection(fun, 0, G.lmax); 
            end
            [ filter_bank, shifted_ends] = mcsfb_design_filter_bank_no_fourier( G, num_bands, band_ends, param);      
        else
            % 6.spectrum_adapted = 0; spacing = log; check_fourier=0;
            for j = 1:num_bands-1
                band_ends(num_bands-j+1) = G.lmax*(1/2)^j;
            end
            [filter_bank, shifted_ends] = mcsfb_design_filter_bank_no_fourier(G, num_bands, band_ends, param);
        end  
    else
        if param.spectrum_adapted
            % 7. spectrum_adapted=1; spacing = even; check_fourier=0
            for j = 1:num_bands-1
                fun = @(x) (G.spectrum_cdf_approx(x) - j/num_bands);
                band_ends(j+1) = bisection(fun, 0, G.lmax); 
            end
            [filter_bank, shifted_ends] = mcsfb_design_filter_bank_no_fourier(G, num_bands, band_ends, param);
        else
            % 8. spectrum_adapted = 0; spacing = even; check_fourier=0;
            for j = 1:num_bands-1
                band_ends(j+1) = G.lmax/num_bands*j;
            end
            [filter_bank, shifted_ends] = mcsfb_design_filter_bank_no_fourier(G, num_bands, band_ends, param);
        end
    end
    
end
end