%% 生成仿真所需的 data_TJSP.mat 文件
% 该文件包含以下变量：
%   GMsettings.Location_Tx : 基站初始位置 [N_bs x 3] (x,y,z)
%   XY                     : 建筑物轮廓点 [N_points x 2] (x,y)
%   building_start_index   : 每个建筑物在 XY 中的起始索引 [N_buildings x 1]
%   building_end_index     : 每个建筑物在 XY 中的结束索引 [N_buildings x 1]
%   num_of_buildings       : 建筑物数量
%   main_walls             : 墙体三角形列表 [N_walls x 12] (4个角点，每个角点 x,y,z)

clear; clc;

%% 1. 定义校园范围（示例：假设校园为 1000m x 800m 的矩形）
x_min = 0;   x_max = 1000;
y_min = 0;   y_max = 800;

%% 2. 生成模拟的基站初始位置（例如 3 个基站）
num_BS = 3;
GMsettings.Location_Tx = [
    200, 200, 10;    % 基站1 (x, y, z)
    500, 400, 10;    % 基站2
    800, 600, 10     % 基站3
];

%% 3. 生成模拟的建筑物轮廓
% 每个建筑物由顺时针或逆时针的多边形顶点构成（闭合，首尾不重复）
% XY 为所有建筑物顶点的集合，building_start_index/end_index 指出每个建筑物的范围

% 建筑物1：矩形 (200,200) 到 (350,300)
buildings{1}.xy = [200,200; 350,200; 350,300; 200,300];
% 建筑物2：矩形 (600,500) 到 (750,600)
buildings{2}.xy = [600,500; 750,500; 750,600; 600,600];
% 建筑物3：L形区域（示例）
buildings{3}.xy = [400,600; 480,600; 480,680; 550,680; 550,750; 400,750];

num_of_buildings = length(buildings);

% 将所有建筑物的顶点顺序存入 XY
XY = [];
building_start_index = zeros(num_of_buildings,1);
building_end_index   = zeros(num_of_buildings,1);

for i = 1:num_of_buildings
    building_start_index(i) = size(XY,1) + 1;
    XY = [XY; buildings{i}.xy];
    building_end_index(i)   = size(XY,1);
end

%% 4. 生成墙体三角形（main_walls）用于精确遮挡检测
% 每个墙体由四个角点（三维）构成，将被分割为两个三角形。
% 假设所有建筑物墙体高度为 5 米（从地面 z=0 到 z=5）。
%
% main_walls 是一个矩阵，每行有 12 个值：
%  [x1 y1 z1, x2 y2 z2, x3 y3 z3, x4 y4 z4]
% 表示一个矩形墙面的四个角点（顺序任意，但需能构成两个三角形）。

wall_height = 5;  % 墙体高度
main_walls = [];

for i = 1:num_of_buildings
    verts = buildings{i}.xy;
    n_verts = size(verts,1);
    for j = 1:n_verts
        % 当前顶点与下一个顶点构成一条棱边
        p1 = verts(j,:);
        p2 = verts(mod(j,n_verts)+1, :);
        
        % 墙体的四个角点（底面两个点，顶面两个点）
        corner1 = [p1, 0];
        corner2 = [p2, 0];
        corner3 = [p2, wall_height];
        corner4 = [p1, wall_height];
        
        % 按顺序存储四个角点 (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        wall = [corner1, corner2, corner3, corner4];
        main_walls = [main_walls; wall];
    end
end

% （可选）添加一些独立的墙体（如围墙）也可以按同样方式添加

%% 5. 保存到 data_TJSP.mat
save_path = 'D:\Propagation_code\Scenario_data\data_TJSP.mat';
% 确保目标文件夹存在
folder = fileparts(save_path);
if ~exist(folder, 'dir')
    mkdir(folder);
 
end

save(save_path, 'GMsettings', 'XY', 'building_start_index', 'building_end_index', ...
                'num_of_buildings', 'main_walls');

fprintf('✅ data_TJSP.mat 已成功生成于：\n%s\n', save_path);
fprintf('包含变量：\n');
fprintf('  GMsettings.Location_Tx (%dx3)\n', size(GMsettings.Location_Tx,1));
fprintf('  XY (%dx2)\n', size(XY,1));
fprintf('  building_start_index (%d)\n', length(building_start_index));
fprintf('  building_end_index (%d)\n', length(building_end_index));
fprintf('  num_of_buildings = %d\n', num_of_buildings);
fprintf('  main_walls (%dx12)\n', size(main_walls,1));