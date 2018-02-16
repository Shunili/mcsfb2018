function [filter_bank, band_ends, shifted_ends] = mcsfb_design_filter_bank(G, num_bands,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(param, 'spacing')
    param.spacing = 'logarithm'; % the other option is even
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
shifted_ends(num_bands+1) = G.lmax;


if param.spacing == 'logarithm'
    if gsp_check_fourier(G) %minima of cdf
        for k = 2:num_bands
            idx = G.N*(1/2)^k;
            shifted_ends(k) = G.e(idx);
        end

        for l = 1:num_bands
            filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
        end 
    else
        if param.spectrum_adapted == 1%minima of cdf % put four cases (if spectrum_adapated=0, don't need to approximate full spectrum)
            if ~isfield(G,'spectrum_cdf_approx')
                G=gsp_spectrum_cdf_approx(G);
            end
            %G.spectrum_cdf_approx = @(x) tern(G, x); %redefine cdf such that cdf = 0 for x<0  

            band_ends = zeros(num_bands+1,1);
            band_ends(1) = 0;
            band_ends(num_bands+1) = G.lmax;

            % compute pdf
            xx = 0:0.001:G.lmax;
            delta=.001;
            G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative

            % approximate eigenvalues
            for j = 1:num_bands-1
                fun = @(x) (G.spectrum_cdf_approx(x)-1/(2^j));
                band_ends(length(band_ends)-j) = bisection(fun, 0, G.lmax); 
            end

            if strcmp(param.band_structure, 'method 1') %minima of cdf
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

            elseif strcmp(param.band_structure, 'method 2') %search in the interval

                % search the minima in range 
                for k = 2:num_bands
                    shifted_ends(k) = fminbnd(G.spectrum_pdf_approx,(band_ends(k)+band_ends(k-1))/2, (band_ends(k+1)+band_ends(k))/2);
                end

                for l = 1:num_bands
                    filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
                end 
            else
                error('Unknown graph type');
            end
        end
    end
end
    
if param.spacing == 'even'
    if gsp_check_fourier(G) %minima of cdf
        for k = 2:num_bands
            idx = G.N*(1/2)^k;
            shifted_ends(k) = G.e(idx);
        end

        for l = 1:num_bands
            filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
        end 
    else
        if param.spectrum_adapted == 1%minima of cdf % put four cases (if spectrum_adapated=0, don't need to approximate full spectrum)
            if ~isfield(G,'spectrum_cdf_approx')
                G=gsp_spectrum_cdf_approx(G);
            end
            %G.spectrum_cdf_approx = @(x) tern(G, x); %redefine cdf such that cdf = 0 for x<0  

            band_ends = zeros(num_bands+1,1);
            band_ends(1) = 0;
            band_ends(num_bands+1) = G.lmax;

            % compute pdf
            xx = 0:0.001:G.lmax;
            delta=.001;
            G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);% first derivative

            % approximate eigenvalues
            for j = 1:num_bands-1
                fun = @(x) (G.spectrum_cdf_approx(x)-1/(2^j));
                band_ends(length(band_ends)-j) = bisection(fun, 0, G.lmax); 
            end

            if strcmp(param.band_structure, 'method 1') %minima of cdf
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

            elseif strcmp(param.band_structure, 'method 2') %search in the interval

                % search the minima in range 
                for k = 2:num_bands
                    shifted_ends(k) = fminbnd(G.spectrum_pdf_approx,(band_ends(k)+band_ends(k-1))/2, (band_ends(k+1)+band_ends(k))/2);
                end

                for l = 1:num_bands
                    filter_bank{l}=@(x) ((shifted_ends(l) <= x) & (x <= shifted_ends(l+1)));
                end 
            else
                error('Unknown graph type');
            end
        end
    end
end

