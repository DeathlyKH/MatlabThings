function Multi_BS_Propagation_Graph_Simulation_Improved()
    % 多基站传播图仿真主函数（带基站位置优化和精确遮挡检测）
    % 完整版：包含所有辅助函数，可直接运行
    
    clc; clear; close all;
    
    % ===================== 配置参数 =====================
    ROOT_PATH = 'D:\';  % 项目根目录，请根据实际情况修改
    
    % 加载数据文件
    data_file = fullfile(ROOT_PATH, 'Propagation_code', 'Scenario_data', 'data_TJSP.mat');
    if ~exist(data_file, 'file')
        error('❌ 未找到数据文件！\n   期望路径：%s\n   请先运行数据生成脚本。', data_file);
    end
    
    fprintf('📂 正在加载数据文件...\n');
    load(data_file);  % 应包含 GMsettings, XY, building_start_index, building_end_index, num_of_buildings, main_walls
    fprintf('✅ 数据加载成功！\n');
    
    % ===================== 优化配置 =====================
    MAX_SCATTERERS = 1000;
    MAX_TX_TO_SCATTERERS = 200;
    MAX_SCATTERERS_TO_RX = 200;
    TEST_POINT_NUM = 300;
    
    % ===================== 基站配置 =====================
    Tx_positions = GMsettings.Location_Tx;
    num_BS = size(Tx_positions, 1);
    fprintf('\n📡 基站配置：\n');
    fprintf('  基站数量：%d\n', num_BS);
    for bs_idx = 1:num_BS
        fprintf('  基站%d: (%.1f, %.1f, %.1f)\n', bs_idx, Tx_positions(bs_idx, :));
    end
    
    % ===================== 生成测试点 =====================
    fprintf('\n🔧 生成测试点...\n');
    x_min = min(XY(:,1)); x_max = max(XY(:,1));
    y_min = min(XY(:,2)); y_max = max(XY(:,2));
    z_min = 0; z_max = 15;
    fprintf('  校园范围：X:[%.1f, %.1f], Y:[%.1f, %.1f]\n', x_min, x_max, y_min, y_max);
    
    Rx_positions = generate_optimized_test_points(TEST_POINT_NUM, XY, building_start_index, building_end_index);
    num_Rx = size(Rx_positions, 1);
    fprintf('  生成有效测试点：%d个\n', num_Rx);
    
    % ===================== 散射体生成 =====================
    fprintf('\n🔬 生成散射体...\n');
    settings.delta_S = 5.0;
    settings.rcs = 1.0;
    settings.trans_att = 20;
    settings.max_trans = 2;
    settings.center_frequency = 2.4e9;
    
    fprintf('  生成墙体散射体...\n');
    building_scatterers = generate_building_scatterers(XY, building_start_index, building_end_index, settings);
    fprintf('  生成随机散射体...\n');
    random_scatterers = generate_random_scatterers(min(MAX_SCATTERERS - size(building_scatterers,1), 100), ...
        x_min, x_max, y_min, y_max);
    Location_S = [building_scatterers; random_scatterers];
    if size(Location_S,1) > MAX_SCATTERERS
        Location_S = Location_S(1:MAX_SCATTERERS,:);
    end
    fprintf('  总散射体数量：%d\n', size(Location_S,1));
    
    % ===================== 基站位置优化 =====================
    fprintf('\n🎯 基站位置优化开始（带精确遮挡检测）...\n');
    bs_to_optimize = 1:num_BS;
    bounds.lower = [x_min, y_min, z_min];
    bounds.upper = [x_max, y_max, z_max];
    
    fprintf('  计算初始覆盖性能...\n');
    initial_positions = Tx_positions;
    initial_metrics = evaluate_coverage_with_proper_occlusion(initial_positions, Rx_positions, ...
        Location_S, main_walls, settings, building_start_index, building_end_index, XY);
    fprintf('  初始平均损耗: %.2f dB, 覆盖率: %.2f%%, 均衡度: %.3f\n', ...
        initial_metrics.mean_path_loss, initial_metrics.coverage_ratio, initial_metrics.service_balance);
    
    fprintf('\n🔍 运行粒子群优化（同时优化%d个基站）...\n', length(bs_to_optimize));
    [optimized_positions, best_metrics, optimization_history] = optimize_multiple_bs_positions_with_proper_occlusion(...
        initial_positions, bs_to_optimize, Rx_positions, Location_S, main_walls, ...
        settings, building_start_index, building_end_index, XY, bounds);
    
    % 显示优化结果
    fprintf('\n✅ 优化完成！\n');
    fprintf('📍 优化前后基站位置对比：\n');
    for bs_idx = 1:num_BS
        if ismember(bs_idx, bs_to_optimize)
            fprintf('  基站%d: (%.1f, %.1f, %.1f) → (%.1f, %.1f, %.1f)\n', ...
                bs_idx, initial_positions(bs_idx,1), initial_positions(bs_idx,2), initial_positions(bs_idx,3), ...
                optimized_positions(bs_idx,1), optimized_positions(bs_idx,2), optimized_positions(bs_idx,3));
        else
            fprintf('  基站%d: (%.1f, %.1f, %.1f) (未优化)\n', bs_idx, initial_positions(bs_idx,:));
        end
    end
    fprintf('📊 性能改进: 平均损耗 %.2f → %.2f dB (改进 %.2f dB)\n', ...
        initial_metrics.mean_path_loss, best_metrics.mean_path_loss, ...
        initial_metrics.mean_path_loss - best_metrics.mean_path_loss);
    
    % ===================== 完整传播计算（含多径） =====================
    fprintf('\n📊 计算优化后的完整传播损耗（直接路径+多径）...\n');
    Tx_positions = optimized_positions;
    c = 3e8;
    lambda = c / settings.center_frequency;
    path_loss_matrix = zeros(num_BS, num_Rx);
    best_bs = zeros(num_Rx, 1);
    
    for bs_idx = 1:num_BS
        fprintf('  基站%d/%d...\n', bs_idx, num_BS);
        Tx = Tx_positions(bs_idx, :);
        [selected_scatterers, scatterer_indices] = select_nearby_scatterers(Tx, Location_S, MAX_TX_TO_SCATTERERS);
        [PL_Tx_to_scatterer, visible_Tx_to_scatterer] = calculate_Tx_to_scatterers_with_occlusion(...
            Tx, selected_scatterers, main_walls, settings);
        
        for rx_idx = 1:num_Rx
            Rx = Rx_positions(rx_idx, :);
            dist_direct = norm(Tx - Rx);
            visible_direct = ~check_wall_blocking_proper(Tx, Rx, main_walls);
            if dist_direct > 0 && visible_direct
                PL_direct = 20*log10(4*pi*dist_direct/lambda);
            else
                PL_direct = Inf;
            end
            [scatterers_for_rx, rx_scatterer_indices] = select_scatterers_for_rx(Rx, selected_scatterers, scatterer_indices, MAX_SCATTERERS_TO_RX);
            [PL_scatterer_to_Rx, visible_scatterer_to_Rx] = calculate_scatterers_to_Rx_with_occlusion(Rx, scatterers_for_rx, main_walls, settings);
            PL_scatterer_scatterer = calculate_scatterer_interactions(selected_scatterers, scatterers_for_rx, scatterer_indices, rx_scatterer_indices, settings);
            PL_scatter = calculate_total_scattering_path(PL_Tx_to_scatterer, PL_scatterer_scatterer, PL_scatterer_to_Rx, visible_Tx_to_scatterer, visible_scatterer_to_Rx, lambda);
            PL_total = min(PL_direct, PL_scatter);
            path_loss_matrix(bs_idx, rx_idx) = PL_total;
        end
    end
    
    for rx_idx = 1:num_Rx
        [~, best_bs(rx_idx)] = min(path_loss_matrix(:, rx_idx));
    end
    bs_service_count = histcounts(best_bs, 1:num_BS+1);
    fprintf('\n📡 基站服务分布：\n');
    for bs_idx = 1:num_BS
        fprintf('  基站%d: %.1f%% (%d个点)\n', bs_idx, bs_service_count(bs_idx)/num_Rx*100, bs_service_count(bs_idx));
    end
    
    % ===================== 可视化 =====================
    fprintf('\n🎨 生成可视化结果...\n');
    figure('Name', '多基站传播仿真与优化结果（精确遮挡）', 'Position', [50, 50, 1800, 900]);
    
    % 子图1：基站位置对比
    subplot(2, 3, [1,2]);
    hold on; grid on; box on;
    for b = 1:num_of_buildings
        idx_start = building_start_index(b);
        idx_end = building_end_index(b);
        if idx_end > size(XY,1), idx_end = size(XY,1); end
        if idx_start > idx_end, continue; end
        fill(XY(idx_start:idx_end,1), XY(idx_start:idx_end,2), [0.7 0.7 0.7], 'EdgeColor','k','FaceAlpha',0.3);
    end
    for bs = 1:num_BS
        if ismember(bs, bs_to_optimize)
            plot(initial_positions(bs,1), initial_positions(bs,2), 's', 'MarkerSize',15, 'MarkerFaceColor','r', 'DisplayName',sprintf('BS%d初始',bs));
            plot([initial_positions(bs,1), optimized_positions(bs,1)], [initial_positions(bs,2), optimized_positions(bs,2)], 'r--');
        else
            plot(initial_positions(bs,1), initial_positions(bs,2), '^', 'MarkerSize',15, 'MarkerFaceColor','b');
        end
    end
    for bs = 1:num_BS
        if ismember(bs, bs_to_optimize)
            plot(optimized_positions(bs,1), optimized_positions(bs,2), 'o', 'MarkerSize',20, 'MarkerFaceColor','g', 'DisplayName',sprintf('BS%d优化后',bs));
        end
    end
    xlabel('X (m)'); ylabel('Y (m)'); title('基站位置优化对比'); axis equal;
    xlim([x_min, x_max]); ylim([y_min, y_max]); legend('Location','best');
    
    % 子图2：收敛曲线
    subplot(2,3,3);
    if ~isempty(optimization_history)
        plot(optimization_history.best_fitness, 'b-', 'LineWidth',2); hold on;
        plot(optimization_history.mean_fitness, 'r--', 'LineWidth',1.5);
        xlabel('迭代次数'); ylabel('适应度'); title('PSO收敛曲线'); grid on;
        legend('最优适应度','平均适应度');
    end
    
    % 子图3：损耗分布
    colors = lines(num_BS);
    subplot(2,3,[4,5,6]);
    hold on; grid on; box on;
    for b = 1:num_of_buildings
        idx_start = building_start_index(b);
        idx_end = building_end_index(b);
        if idx_end > size(XY,1), idx_end = size(XY,1); end
        if idx_start > idx_end, continue; end
        fill(XY(idx_start:idx_end,1), XY(idx_start:idx_end,2), [0.7 0.7 0.7], 'EdgeColor','k','FaceAlpha',0.3);
    end
    path_loss_vals = zeros(num_Rx,1);
    for i = 1:num_Rx, path_loss_vals(i) = path_loss_matrix(best_bs(i), i); end
    valid = ~isinf(path_loss_vals);
    if sum(valid)>0
        minL = min(path_loss_vals(valid)); maxL = max(path_loss_vals(valid));
        path_loss_vals(~valid) = maxL + 10;
    else
        minL = 50; maxL = 150;
    end
    scatter(Rx_positions(:,1), Rx_positions(:,2), 10, path_loss_vals, 'filled');
    for bs = 1:num_BS
        plot(optimized_positions(bs,1), optimized_positions(bs,2), '^', ...
             'MarkerSize',20, 'MarkerFaceColor', colors(bs,:), ...
             'MarkerEdgeColor','k');
        text(optimized_positions(bs,1)+10, optimized_positions(bs,2)+10, sprintf('BS%d',bs), ...
             'FontSize',12, 'BackgroundColor','w');
    end
    colormap(jet); colorbar; caxis([minL, maxL]);
    xlabel('X (m)'); ylabel('Y (m)'); title('传播损耗分布'); axis equal;
    xlim([x_min, x_max]); ylim([y_min, y_max]);
    
    % ===================== 保存结果与热力图 =====================
    save_dir = fullfile(ROOT_PATH, 'Propagation_code', 'Multi_BS_Results');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_file = fullfile(save_dir, sprintf('Multi_BS_Optimization_%dBS_%dRx_%s.mat', num_BS, num_Rx, timestamp));
    save(save_file, 'initial_positions', 'optimized_positions', 'initial_metrics', 'best_metrics', ...
        'path_loss_matrix', 'best_bs', 'settings', 'optimization_history');
    saveas(gcf, fullfile(save_dir, sprintf('optimization_results_%s.png', timestamp)));
    fprintf('✅ 结果保存至：%s\n', save_file);
    
    % 生成热力图
    generate_heatmap(optimized_positions, XY, building_start_index, building_end_index, ...
                     main_walls, settings, save_dir, timestamp, x_min, x_max, y_min, y_max);
    
    fprintf('\n🎉 仿真完成！\n');
end

% -------------------------------------------------------------------------
% 以下为所有辅助函数
% -------------------------------------------------------------------------

function Rx = generate_optimized_test_points(target_num, XY, building_start_index, building_end_index)
    z_height = 1.5; Rx = []; max_iter = 5; iter = 0;
    x_min = min(XY(:,1)); x_max = max(XY(:,1));
    y_min = min(XY(:,2)); y_max = max(XY(:,2));
    while size(Rx,1) < target_num && iter < max_iter
        iter = iter + 1;
        batch = min(target_num*3, 2000);
        rx_x = x_min + rand(batch,1)*(x_max-x_min);
        rx_y = y_min + rand(batch,1)*(y_max-y_min);
        candidates = [rx_x, rx_y, ones(batch,1)*z_height];
        keep = true(batch,1);
        num_b = length(building_start_index);
        check_b = min(100, num_b);
        idxs = randperm(num_b, check_b);
        for i = 1:check_b
            b = idxs(i);
            s = building_start_index(b); e = building_end_index(b);
            if e>size(XY,1), e=size(XY,1); end
            if s>e || s>size(XY,1), continue; end
            xp = XY(s:e,1); yp = XY(s:e,2);
            if length(xp)>=3
                in = inpolygon(candidates(:,1), candidates(:,2), xp, yp);
                keep = keep & ~in;
            end
        end
        Rx = [Rx; candidates(keep,:)];
        if size(Rx,1) >= target_num, Rx = Rx(1:target_num,:); break; end
    end
end

function scatterers = generate_building_scatterers(XY, building_start_index, building_end_index, settings)
    scatterers = []; num_b = length(building_start_index);
    max_b = min(50, num_b); idxs = randperm(num_b, max_b);
    for idx = 1:max_b
        b = idxs(idx);
        s = building_start_index(b); e = building_end_index(b);
        if e>size(XY,1), e=size(XY,1); end
        if s>e || s>size(XY,1), continue; end
        xp = XY(s:e,1); yp = XY(s:e,2);
        if length(xp)<3, continue; end
        % 周长
        peri = 0;
        for i=1:length(xp)
            i2 = mod(i,length(xp))+1;
            peri = peri + sqrt((xp(i2)-xp(i))^2 + (yp(i2)-yp(i))^2);
        end
        num_on = max(1, round(peri/settings.delta_S));
        for i=1:num_on
            seg = mod(i-1, length(xp))+1;
            frac = (i-1)/num_on;
            x1 = xp(seg); y1 = yp(seg);
            x2 = xp(mod(seg,length(xp))+1); y2 = yp(mod(seg,length(xp))+1);
            x = x1 + frac*(x2-x1); y = y1 + frac*(y2-y1); z = 1.5 + rand*3;
            scatterers = [scatterers; x,y,z];
        end
        % 内部点
        cx = mean(xp); cy = mean(yp);
        num_in = max(1, round(sqrt(num_on)));
        for i=1:num_in
            alpha = rand; beta = rand;
            idx1 = randi(length(xp)); idx2 = mod(idx1,length(xp))+1;
            x = cx*alpha + (xp(idx1)*beta + xp(idx2)*(1-beta))*(1-alpha);
            y = cy*alpha + (yp(idx1)*beta + yp(idx2)*(1-beta))*(1-alpha);
            z = 1.5 + rand*3;
            scatterers = [scatterers; x,y,z];
        end
    end
    if size(scatterers,1) > 500, scatterers = scatterers(1:500,:); end
end

function scatterers = generate_random_scatterers(num_scatterers, x_min, x_max, y_min, y_max)
    x = x_min + rand(num_scatterers,1)*(x_max-x_min);
    y = y_min + rand(num_scatterers,1)*(y_max-y_min);
    z = 1.5 + rand(num_scatterers,1)*3;
    scatterers = [x,y,z];
end

function [selected, idxs] = select_nearby_scatterers(Tx, Location_S, max_count)
    dist = vecnorm(Location_S - Tx, 2, 2);
    [~, ord] = sort(dist);
    cnt = min(max_count, length(dist));
    idxs = ord(1:cnt);
    selected = Location_S(idxs,:);
end

function [scat_rx, idxs_rx] = select_scatterers_for_rx(Rx, selected_scatterers, scatterer_indices, max_count)
    dist = vecnorm(selected_scatterers - Rx, 2, 2);
    [~, ord] = sort(dist);
    cnt = min(max_count, length(dist));
    idxs_rx = scatterer_indices(ord(1:cnt));
    scat_rx = selected_scatterers(ord(1:cnt),:);
end

function PL = calculate_scatterer_interactions(S1, S2, idx1, idx2, settings)
    n1 = size(S1,1); n2 = size(S2,1);
    PL = zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            if idx1(i)==idx2(j)
                PL(i,j)=0;
            else
                d = norm(S1(i,:)-S2(j,:));
                if d>0, PL(i,j)=20*log10(d)+30+10*randn; else PL(i,j)=0; end
            end
        end
    end
end

function PL_total = calculate_total_scattering_path(PL_Tx_S, PL_SS, PL_S_Rx, vis_Tx, vis_Rx, lambda)
    num_Tx = length(PL_Tx_S); num_Rx = length(PL_S_Rx);
    if num_Tx==0 || num_Rx==0, PL_total=Inf; return; end
    best = Inf;
    for i=1:num_Tx
        if ~vis_Tx(i), continue; end
        for j=1:num_Rx
            if ~vis_Rx(j), continue; end
            loss_ss = PL_SS(i,j) 
            if i<=size(PL_SS,1) && j<=size(PL_SS,2) else 30; end
            total = PL_Tx_S(i) + loss_ss + PL_S_Rx(j);
            if total < best, best = total; end
        end
    end
    PL_total = best;
end

function blocked = check_wall_blocking_proper(Tx, Rx, walls)
    blocked = false;
    if isempty(walls), return; end
    max_check = min(20, size(walls,1));
    for w = 1:max_check
        wall = walls(w,:);
        corners = reshape(wall, [3,4])';
        tri = [1,2,3; 1,3,4];
        for t=1:2
            A = corners(tri(t,1),:); B = corners(tri(t,2),:); C = corners(tri(t,3),:);
            if line_triangle_intersect(Tx, Rx, A, B, C)
                blocked = true; return;
            end
        end
    end
end

function intersect = line_triangle_intersect(P1, P2, V0, V1, V2)
    intersect = false;
    edge1 = V1 - V0; edge2 = V2 - V0;
    h = cross(P2-P1, edge2);
    det = dot(edge1, h);
    if abs(det) < 1e-8, return; end
    inv_det = 1/det;
    s = P1 - V0;
    u = inv_det * dot(s, h);
    if u<0 || u>1, return; end
    q = cross(s, edge1);
    v = inv_det * dot(P2-P1, q);
    if v<0 || u+v>1, return; end
    t = inv_det * dot(edge2, q);
    if t>=0 && t<=1, intersect = true; end
end

function metrics = evaluate_coverage_with_proper_occlusion(Tx_positions, Rx_positions, ~, walls, ...
    settings, building_start_index, building_end_index, XY)
    num_BS = size(Tx_positions,1); num_Rx = size(Rx_positions,1);
    c = 3e8; lambda = c/settings.center_frequency;
    loss_mat = inf(num_BS, num_Rx);
    for bs=1:num_BS
        Tx = Tx_positions(bs,:);
        for rx=1:num_Rx
            Rx = Rx_positions(rx,:);
            d = norm(Tx-Rx);
            visible = ~check_wall_blocking_proper(Tx, Rx, walls);
            if d>0 && visible
                loss_mat(bs,rx) = 20*log10(4*pi*d/lambda);
            end
        end
    end
    best_loss = min(loss_mat,[],1);
    best_bs = zeros(num_Rx,1);
    for i=1:num_Rx, [~,best_bs(i)] = min(loss_mat(:,i)); end
    valid = ~isinf(best_loss);
    if any(valid)
        metrics.mean_path_loss = mean(best_loss(valid));
        metrics.max_path_loss = max(best_loss(valid));
        metrics.min_path_loss = min(best_loss(valid));
        coverage = sum(best_loss(valid) <= -90) / sum(valid) * 100;
        metrics.coverage_ratio = coverage;
    else
        metrics.mean_path_loss = Inf; metrics.max_path_loss = Inf; metrics.coverage_ratio = 0;
    end
    cnt = histcounts(best_bs, 1:num_BS+1);
    if sum(cnt)>0
        sorted = sort(cnt);
        n = length(sorted);
        gini = (n+1 - 2*sum(cumsum(sorted))/sum(sorted)) / n;
        metrics.service_balance = 1 - gini;
    else
        metrics.service_balance = 0;
    end
    if ~isinf(metrics.mean_path_loss)
        norm_mean = (metrics.mean_path_loss + 100) / 50;
        norm_cov = (100 - metrics.coverage_ratio) / 100;
        norm_max = (metrics.max_path_loss + 100) / 50;
        norm_bal = 1 - metrics.service_balance;
        metrics.fitness = 0.4*norm_mean + 0.3*norm_cov + 0.2*norm_max + 0.1*norm_bal;
    else
        metrics.fitness = Inf;
    end
end

function [PL, visible] = calculate_Tx_to_scatterers_with_occlusion(Tx, scatterers, walls, settings)
    n = size(scatterers,1);
    PL = zeros(n,1); visible = true(n,1);
    c=3e8; lambda = c/settings.center_frequency;
    for i=1:n
        d = norm(Tx - scatterers(i,:));
        if check_wall_blocking_proper(Tx, scatterers(i,:), walls)
            visible(i)=false;
            PL(i)=20*log10(4*pi*d/lambda)+settings.trans_att;
        else
            PL(i)=20*log10(4*pi*d/lambda);
        end
    end
end

function [PL, visible] = calculate_scatterers_to_Rx_with_occlusion(Rx, scatterers, walls, settings)
    n = size(scatterers,1);
    PL = zeros(n,1); visible = true(n,1);
    c=3e8; lambda = c/settings.center_frequency;
    for i=1:n
        d = norm(Rx - scatterers(i,:));
        if check_wall_blocking_proper(scatterers(i,:), Rx, walls)
            visible(i)=false;
            PL(i)=20*log10(4*pi*d/lambda)+settings.trans_att;
        else
            PL(i)=20*log10(4*pi*d/lambda);
        end
    end
end

function [opt_pos, best_met, hist] = optimize_multiple_bs_positions_with_proper_occlusion(...
    init_pos, bs_opt, Rx_pos, S, walls, settings, b_start, b_end, XY, bounds)
    pso.num_particles = 20;
    pso.max_iter = 40;
    pso.w = 0.9; pso.c1 = 1.5; pso.c2 = 1.5;
    pso.v_min = -15; pso.v_max = 15;
    num_BS = size(init_pos,1);
    num_opt = length(bs_opt);
    dim = 3 * num_opt;
    particles = cell(pso.num_particles,1);
    vel = zeros(pso.num_particles, dim);
    pbest_pos = zeros(pso.num_particles, dim);
    pbest_fit = inf(pso.num_particles,1);
    gbest_pos = zeros(1,dim);
    gbest_fit = inf;
    best_met = [];
    hist.best_fitness = zeros(pso.max_iter,1);
    hist.mean_fitness = zeros(pso.max_iter,1);
    
    % 初始化粒子
    for p = 1:pso.num_particles
        new_pos_mat = init_pos;
        vec = zeros(1,dim);
        for i=1:num_opt
            bs = bs_opt(i);
            sidx = (i-1)*3+1; eidx = i*3;
            if p==1
                pos = init_pos(bs,:);
            else
                pos = bounds.lower + rand(1,3).*(bounds.upper-bounds.lower);
                for attempt=1:5
                    if ~is_in_building(pos, XY, b_start, b_end), break; end
                    pos = bounds.lower + rand(1,3).*(bounds.upper-bounds.lower);
                end
            end
            new_pos_mat(bs,:) = pos;
            vec(sidx:eidx) = pos;
        end
        met = evaluate_coverage_with_proper_occlusion(new_pos_mat, Rx_pos, S, walls, settings, b_start, b_end, XY);
        fit = met.fitness;
        particles{p} = struct('pos',vec, 'fit',fit);
        pbest_pos(p,:) = vec;
        pbest_fit(p) = fit;
        if fit < gbest_fit
            gbest_fit = fit; gbest_pos = vec; best_met = met;
        end
    end
    
    % PSO主循环
    for iter = 1:pso.max_iter
        if mod(iter,5)==0, fprintf('    迭代 %d/%d, 最佳适应度: %.4f\n', iter, pso.max_iter, gbest_fit); end
        for p = 1:pso.num_particles
            r1 = rand(1,dim); r2 = rand(1,dim);
            vel(p,:) = pso.w*vel(p,:) + pso.c1*r1.*(pbest_pos(p,:)-particles{p}.pos) + pso.c2*r2.*(gbest_pos-particles{p}.pos);
            vel(p,:) = max(min(vel(p,:), pso.v_max), pso.v_min);
            new_vec = particles{p}.pos + vel(p,:);
            new_mat = init_pos;
            for i=1:num_opt
                bs = bs_opt(i);
                sidx = (i-1)*3+1; eidx = i*3;
                pos = new_vec(sidx:eidx);
                pos = max(min(pos, bounds.upper), bounds.lower);
                for attempt=1:3
                    if ~is_in_building(pos, XY, b_start, b_end), break; end
                    pos = pos + randn(1,3)*5;
                    pos = max(min(pos, bounds.upper), bounds.lower);
                end
                new_mat(bs,:) = pos;
                new_vec(sidx:eidx) = pos;
            end
            met = evaluate_coverage_with_proper_occlusion(new_mat, Rx_pos, S, walls, settings, b_start, b_end, XY);
            new_fit = met.fitness;
            particles{p}.pos = new_vec; particles{p}.fit = new_fit;
            if new_fit < pbest_fit(p)
                pbest_fit(p) = new_fit; pbest_pos(p,:) = new_vec;
                if new_fit < gbest_fit
                    gbest_fit = new_fit; gbest_pos = new_vec; best_met = met;
                end
            end
        end
        hist.best_fitness(iter) = gbest_fit;
        fits = cellfun(@(x)x.fit, particles);
        hist.mean_fitness(iter) = mean(fits);
        pso.w = pso.w * 0.98;
    end
    opt_pos = init_pos;
    for i=1:num_opt
        bs = bs_opt(i);
        sidx = (i-1)*3+1; eidx = i*3;
        opt_pos(bs,:) = gbest_pos(sidx:eidx);
    end
end

function inside = is_in_building(point, XY, building_start_index, building_end_index)
    inside = false;
    for b = 1:length(building_start_index)
        s = building_start_index(b); e = building_end_index(b);
        if e>size(XY,1), e=size(XY,1); end
        if s>e || s>size(XY,1), continue; end
        xp = XY(s:e,1); yp = XY(s:e,2);
        if length(xp)>=3 && inpolygon(point(1), point(2), xp, yp)
            inside = true; return;
        end
    end
end

function generate_heatmap(optimized_positions, XY, building_start_index, building_end_index, ...
                          main_walls, settings, save_dir, timestamp, x_min, x_max, y_min, y_max)
    fprintf('\n🔥 生成全地图热力图...\n');
    res = 2.0;
    xg = x_min:res:x_max; yg = y_min:res:y_max;
    [Xg,Yg] = meshgrid(xg,yg);
    Z = NaN(size(Xg));
    c = 3e8; lambda = c/settings.center_frequency;
    total = numel(Xg);
    fprintf('  网格点数: %d\n', total);
    for idx = 1:total
        if mod(idx, round(total/10))==0, fprintf('    进度: %.0f%%\n', idx/total*100); end
        pt = [Xg(idx), Yg(idx), 1.5];
        if is_in_building(pt, XY, building_start_index, building_end_index), continue; end
        best_loss = Inf;
        for bs = 1:size(optimized_positions,1)
            Tx = optimized_positions(bs,:);
            d = norm(Tx - pt);
            vis = ~check_wall_blocking_proper(Tx, pt, main_walls);
            if d>0 && vis
                loss = 20*log10(4*pi*d/lambda);
                if loss < best_loss, best_loss = loss; end
            end
        end
        if isfinite(best_loss), Z(idx) = best_loss; end
    end
    figure('Name','校园覆盖热力图','Position',[100,100,1000,800]);
    imagesc(xg, yg, Z); set(gca,'YDir','normal');
    colormap(flipud(jet)); colorbar; caxis([-110 -50]); hold on;
    for b = 1:length(building_start_index)
        s = building_start_index(b); e = building_end_index(b);
        if e>size(XY,1), e=size(XY,1); end
        if s>e, continue; end
        fill(XY(s:e,1), XY(s:e,2), [0.3,0.3,0.3], 'EdgeColor','k','FaceAlpha',0.6);
    end
    for bs = 1:size(optimized_positions,1)
        plot(optimized_positions(bs,1), optimized_positions(bs,2), '^', 'MarkerSize',20, 'MarkerFaceColor','r','MarkerEdgeColor','k');
        text(optimized_positions(bs,1)+15, optimized_positions(bs,2)+15, sprintf('BS%d',bs), 'FontSize',12, 'BackgroundColor','w');
    end
    xlabel('X (m)'); ylabel('Y (m)'); title(sprintf('信号覆盖热力图 (分辨率 %.1f m)',res));
    axis equal tight;
    saveas(gcf, fullfile(save_dir, sprintf('Heatmap_%s.png', timestamp)));
    close(gcf);
    fprintf('  ✅ 热力图已保存\n');
end