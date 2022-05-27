%% LINE TRACKING

x = linspace(-1, 1, 41); y = 1.2*ones(size(x));
test2 = [x; y]; t = size(x, 2);

for m = 1:t
    u1 = test2(:, m);
    for j = 1:100
        dist(j) = norm(u1-w_g(:, :, j));
    end
    
    [~, win_val] = min(dist);
    win = [lx(win_val), ly(win_val)];
    d = repmat(win', 1, 100)-lattice;
    H_g = exp(-(sum(d.^2))/(2*(sig_f^2)));

    s = sum(H_g); s2 = 0;
    
    for k = 1:100
        s1 = H_g(k)*(th_g(:,:,k)+A_g(:,:,k)*(u1-w_g(:,:,k)));
        s2 = s2 + s1;
    end      
    
    theta = s2/s;
    x_o = l1*cos(theta(1)) + l2*cos(theta(1) + theta(2));
    y_o = l1*sin(theta(1)) + l2*sin(theta(1) + theta(2));
    v_o = [x_o; y_o];
    th(:, m) = theta;
end

for i = 1:t
    x_Position(i, :) = [0 l1*cos(th(1, i)) l1*cos(th(1, i)) + l2*cos(th(2, i) + th(1, i))];
    y_Position(i, :) = [0 l1*sin(th(1, i)) l1*sin(th(1, i)) + l2*sin(th(2, i) + th(1, i))];
end


figure; plot(test2(1, :), test2(2, :), '-ok'); hold on;

for i = 1:t
    h = plot(x_Position(i, :), y_Position(i, :), 'k');
    pause(0.5);
end

axis equal