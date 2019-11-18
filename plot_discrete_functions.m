function plot_discrete_functions(X,Gf,Gh,F,H)
%     plot_discrete_functions(X,Gf,Gh)
%     Plots the 2d-discrete representation of an instance of PEP for NoLips.
%     Shows the points x_i, as well as the (sub)gradients of f and h
%     
%     Arguments
%         - X : array of size (n_pts, dim) storing the iterates x_i
%         - Gf: array of size (n_pts, dim) storing the gradients of f at the x_i's
%         - Gh: array of size (n_pts, dim) storing the gradients of h at the x_i's
%         - F : array of size (n_pts) storing the values of f
%         - H : array of size (n_pts) storing the values of h
% 

    n_pts = size(X,2);
    dim = size(X,1);
    
    if ~((dim == 1) || (dim == 2))
        error("plot_discrete_functions accepts only dimension 1 or 2");
    end
    
    labels_x = {'x*'};

    for i = 1:n_pts-1
        labels_x{i+1} = sprintf('x%d', i-1);
    end

    if dim == 1

        f_vals = [0; double(F)]';
        plot(X, f_vals, "o");

        hold on
        quiver(X, f_vals, ones(1,n_pts), Gf, 'AutoScale', 'on', 'Color',...
        'blue', 'LineWidth', 1, 'DisplayName', '\nabla f(x_i)');


        quiver(X, f_vals, ones(1,n_pts), Gf, 'AutoScale', 'on', ...
            'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', '\nabla h(x_i)');

  
        text(X, f_vals,labels_x,'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');

        hold off
        legend

    elseif dim == 2


    plot(X(1,:), X(2,:), "o");
    text(X(1,:), X(2,:),labels_x,'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');

    grid on
    axis equal

    hold on

    quiver(X(1,:), X(2,:), Gf(1,:), Gf(2,:), 'AutoScale', 'off', ...
        'LineWidth', 2, 'DisplayName', '\nabla f');

    quiver(X(1,:), X(2,:), Gh(1,:), Gh(2,:), 'AutoScale', 'off', ...
        'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', '\nabla h');



    hold off;
    legend
end

