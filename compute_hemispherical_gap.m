function [T_hemi_gauss, T_hemi_trapz] = compute_hemispherical_gap(data)
% COMPUTE_HEMISPHERICAL_GAP 计算冠层上半球漫射辐射透过率（半球间隙率）
%
% 输入参数:
%   data - 观测数据矩阵。支持两种结构格式：
%          格式1 (N x 3): [观测天顶角(VZA), 观测方位角(VAA), 间隙率(P_gap)]
%          格式2 (N x 2): [观测天顶角(VZA), 方位角平均间隙率(P_gap)]
%          注：角度单位均为度(Degree)。
%
% 输出参数:
%   T_hemi_gauss - 基于5点高斯-勒让德积分求得的半球间隙率 (精度更高，推荐作为首选物理量)
%   T_hemi_trapz - 基于梯形数值积分求得的半球间隙率 (可作为参考对比)
%
% 调用示例:
%   [gauss_res, trapz_res] = compute_hemispherical_gap(my_data_matrix);

    % 1. 数据解析与预处理 (自适应维度)
    if size(data, 2) == 3
        % 格式1: 提取唯一的天顶角并计算方位角平均
        vza_unique = unique(data(:,1));
        P_gap_avg = zeros(length(vza_unique), 1);
        for i = 1:length(vza_unique)
            idx = (data(:,1) == vza_unique(i));
            P_gap_avg(i) = mean(data(idx, 3)); 
        end
    elseif size(data, 2) == 2
        % 格式2: 直接提取并确保按天顶角升序排列
        vza_unique = data(:, 1);
        P_gap_avg = data(:, 2);
        [vza_unique, sort_idx] = sort(vza_unique);
        P_gap_avg = P_gap_avg(sort_idx);
    else
        error('植被物理模型输入错误：请提供 N x 2 或 N x 3 的数据矩阵。');
    end

    % 2. 物理边界条件约束：水平方向(90度)光程趋于无穷，间隙率设为0
    if max(vza_unique) < 90
        vza_full = [vza_unique; 90];
        P_gap_full = [P_gap_avg; 0];
    else
        vza_full = vza_unique;
        P_gap_full = P_gap_avg;
    end

    % 角度转弧度，用于三角函数计算
    theta_rad = deg2rad(vza_full);

    % -------------------------------------------------------------------------
    % 核心算法 1：5点高斯-勒让德积分 (5-point Gauss-Legendre Quadrature)
    % -------------------------------------------------------------------------
    % 标准区间 [-1, 1] 下的高斯节点与权重 (LAI-2000等仪器的底层数学基础)
    x_gauss = [0; 0.5384693101; -0.5384693101; 0.9061798459; -0.9061798459];
    w_gauss = [0.5688888889; 0.4786286705; 0.4786286705; 0.2369268850; 0.2369268850];

    % 积分区间映射：从标准域 [-1, 1] 映射到真实半球空间 [0, pi/2]
    a = 0; 
    b = pi/2;
    t_nodes = ((b - a)/2) .* x_gauss + (b + a)/2; % 映射后的角度节点(弧度)
    W_weights = ((b - a)/2) .* w_gauss;           % 映射后的权重

    % 利用三次样条(Spline)将离散间隙率投影到连续的高斯节点上
    P_gap_gauss_nodes = spline(theta_rad, P_gap_full, t_nodes);
    
    % 物理合理性截断：间隙率必须介于 [0, 1] 之间
    P_gap_gauss_nodes = max(0, min(1, P_gap_gauss_nodes));

    % 在高斯节点上计算冠层透过率的权重通量：2 * P(theta) * sin(theta) * cos(theta)
    f_gauss_nodes = 2 .* P_gap_gauss_nodes .* sin(t_nodes) .* cos(t_nodes);
    
    % 累加得出高斯积分结果
    T_hemi_gauss = sum(W_weights .* f_gauss_nodes);

    % -------------------------------------------------------------------------
    % 核心算法 2：梯形数值积分法 (Trapezoidal Integration)
    % -------------------------------------------------------------------------
    integrand = 2 .* P_gap_full .* sin(theta_rad) .* cos(theta_rad);
    T_hemi_trapz = trapz(theta_rad, integrand);

end