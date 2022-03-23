##Тартыков Лев ИУ7-64Б, 2022 г
pkg load statistics;
echo off all;

function main()
    X =[-2.79,-3.01,-4.07,-2.85, -2.43,-3.20,-3.72,-4.27,-5.48,-2.38,-4.69,-4.34,-5.08,-5.01,-4.08,-4.20,-4.74,-1.88,-3.25,-2.78,-3.56,-3.54,-3.79,-3.18,-5.08,-4.30,-2.86,-2.45,-3.08,-3.22,-2.76,-3.20,-3.33,-4.91,-4.06,-3.81,-3.96,-3.65,-3.77,-4.60,-5.21,-2.67,-1.95,-2.43,-1.73,-2.50,-3.96,-3.75,-2.70,-4.26,-3.42,-4.07,-4.74,-3.00,-4.37,-5.42,-5.00,-4.08,-2.46,-4.33,-4.08,-3.72,-4.09,-2.96,-3.71,-1.51,-3.70,-6.48,-4.26,-4.39,-3.16,-4.63,-2.66,-2.22,-4.79,-2.46,-3.69,-3.35,-2.32,-4.17,-3.85,-4.93,-2.05,-3.15,-3.49,-5.70,-2.53,-3.85,-4.32,-3.37,-3.98,-3.74,-5.28,-2.56,-3.21,-3.10,-3.78,-3.36,-3.32,-2.59,-2.45,-3.34,-3.20,-4.14,-4.00,-4.79,-4.02,-4.58,-4.45,-3.69,-4.53,-3.98,-4.51,-4.44,-3.78,-4.24,-4.00,-2.46,-2.58,-4.04];
    [bins, counts, count_X, delta, Xn, Y_normpdf, Y_normcdf, Y_ecdf, M_min, M_max, X_without_double] = perform_params(X);
    plot_graphs(bins, counts, count_X, delta, Xn, Y_normpdf, Y_normcdf, Y_ecdf, M_min, M_max, X_without_double);
endfunction

function [bins, counts, count_X, delta, Xn, Y_normpdf, Y_normcdf, Y_ecdf, M_min, M_max, X_without_double] = perform_params(X)
    count_X = length(X);
    M_max = max(X)
    M_min = min(X)
    
    R = M_max - M_min;
    MX = find_MX(X, count_X);
    DX = find_DX(X, MX, count_X);
    
    m = find_m(count_X);
    [counts, bins] = hist(X, m);

    delta = R / m;
    sigma = sqrt(DX);
    abs_MX = abs(MX);
    M_min = M_min - abs_MX;
    M_max = M_max + abs_MX;
    Xn = M_min:delta/20:M_max;

    Y_normpdf = density_ndist(Xn, MX, sigma);
    #Y_normpdf = normpdf(Xn, MX, sigma);
    Y_normcdf = form_normcdf(Xn, MX, sigma);
    #Y_normcdf = normcdf(Xn, MX, sigma);
    [count_elem, X_without_double] = count_number_elems(X, count_X, M_min, M_max);
    Y_ecdf = form_y_ecdf(count_elem, count_X, M_min, M_max);
    #Y_ecdf = empirical_cdf(Xn, X);
endfunction

function output_centers(bins, counts)
    for i = 1 : length(counts)
        fprintf("центр интервала: %f; количество значений: %d\n", bins(i), counts(i));
    endfor;
    fprintf("\n");
endfunction

function [MX] = find_MX(X, count_X)
    MX = sum(X) / count_X;
endfunction

function [DX] = find_DX(X, MX, count_X)
    DX = sum((X - MX).^2) / (count_X - 1);
endfunction

function [m] = find_m(count_X)
    m = floor(log2(count_X)) + 2;
endfunction

function [Y_normpdf] = density_ndist(Xn, MX, sigma)
    Y_normpdf = [];
    count_X = length(Xn);
    for i = 1: count_X
        Y_normpdf(i) = 1 / (sqrt(2 * pi) * sigma) * exp(-(Xn(i) - MX).^2 / (2 * sigma.^2));
    endfor
endfunction

function [count_elem, X_graph] = count_number_elems(X, count_X, M_min, M_max)
    X_sort = sort(X);
    count_elem = [];
    X_without_double = [];
    i = 1;
    index_count = 1;
    while (i < count_X)
        is_all = 0;
        temp_value = X_sort(i);
        j = 1;
        while (is_all == 0)
            if (X_sort(i + j) == temp_value)
                j++;
            else
                is_all = 1;
            endif
        endwhile
        X_without_double(index_count) = temp_value;
        count_elem(index_count) = j;
        index_count += 1;
        i += j; 
    endwhile

    if (X_sort(count_X) != X_sort(count_X - 1))
        count_elem(index_count) = 1;
        X_without_double(index_count) = X_sort(count_X);
    endif
    
    X_graph = [];
    X_graph(1) = X_without_double(1) - abs(M_min);
    X_graph(2) = X_without_double(1);
    X_graph(3) = X_graph(2);
    j = 4;
    for i = 2: length(X_without_double)
        X_graph(j) = X_without_double(i);
        X_graph(j + 1) = X_without_double(i);
        j += 2;
    endfor
    M_max
    X_graph(j) = abs(M_max);  
endfunction

function [Y_normcdf] = form_normcdf(Xn, MX, sigma)
    Y_normpdf = @(Xn) ((1./(sqrt(2.*pi).*sigma)).*(exp((-0.5.*(Xn - MX).^2)./(sigma.^2))));
    for k = 1: length(Xn)
        Y_normcdf(k) = integral(Y_normpdf, -inf, Xn(k));
    endfor
endfunction

function [Y_ecdf] = form_y_ecdf(count_elem, count_X)
    Y_ecdf = [];
    len_celem = length(count_elem);

    Y_ecdf(1) = 0;
    Y_ecdf(2) = 0;
    Y_ecdf(3) = count_elem(1) / count_X;
    Y_ecdf(4) = Y_ecdf(3);
    j = 5;
    for i = 2: len_celem
        Y_ecdf(j) = Y_ecdf(j - 1) + count_elem(i) / count_X;
        Y_ecdf(j + 1) = Y_ecdf(j);
        j += 2;
    endfor
endfunction

function plot_graphs(bins, counts, count_X, delta, Xn, Y_normpdf, Y_normcdf, Y_ecdf, M_min, M_max, X_without_double)
    figure;
    subplot(1, 2, 1);
    bar(bins, counts / (count_X * delta), "histc", 'FaceColor', 'blue');
    hold on;
    plot(Xn, Y_normpdf, 'LineWidth', 3, 'Color', 'green');
    xlim([M_min, M_max]);
    hold on;
    
    subplot(1, 2, 2);
    plot(Xn, Y_normcdf, 'LineWidth', 1, 'Color', 'red');
    xlim([M_min, M_max]);
    hold on;
    plot(X_without_double, Y_ecdf, 'LineWidth', 1, 'Color', 'blue');
    #plot(Xn, Y_ecdf, 'LineWidth', 1, 'Color', 'blue');
    xlim([M_min, M_max]);

endfunction

main()