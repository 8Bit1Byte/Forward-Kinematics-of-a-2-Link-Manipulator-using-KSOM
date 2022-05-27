%% CIRCLE TRACKING
t = 0:pi/30:2*pi; test1 = [1.5*cos(t); 1.5*sin(t)];
% Link lenghts
l1 = 1; l2 = 1;

A_g = 0.1*rand(2, 2, 100);
w_g = 0.1*rand(2, 1, 100);
th_g = 0.1*rand(2, 1, 100);

% 2d lattice formation of size 10x10
[lx, ly] = ind2sub([10, 10], 1:100);
lattice = [lx; ly]; iterations = 6000;

for m = 1:length(t)
    u2 = test1(:, m);
    for j = 1:100
        dist(j) = norm(u2-w_g(:, :, j));
    end
    
    [~, win_val] = min(dist);
    win = [lx(win_val), ly(win_val)];
    d = repmat(win', 1, 100)-lattice;
    H_g = exp(-(sum(d.^2))/(2*(sig_f^2)));

    % Corse Action
    s = sum(H_g); s2 = 0;
    
    for k = 1:100
        s1 = H_g(k)*(th_g(:,:,k)+A_g(:,:,k)*(u2-w_g(:,:,k)));
        s2 = s2 + s1;
    end      
    
    theta = s2/s;
    x_o = l1*cos(theta(1)) + l2*cos(theta(1) + theta(2));
    y_o = l1*sin(theta(1)) + l2*sin(theta(1) + theta(2));
    v_o = [x_o; y_o];
    th(:, m) = theta;
end

for i = 1:length(t)
    x_Position(i, :) = [0 l1*cos(th(1, i)) l1*cos(th(1, i)) + l2*cos(th(2, i) + th(1, i))];
    y_Position(i, :) = [0 l1*sin(th(1, i)) l1*sin(th(1, i)) + l2*sin(th(2, i) + th(1, i))];
end


figure; plot(test1(1, :), test1(2, :), '-ok'); hold on; grid();

for i = 1:length(t)
    h = plot(x_Position(i, :), y_Position(i, :), 'k');
    axis equal
    pause(2.1);
end

axis equal
