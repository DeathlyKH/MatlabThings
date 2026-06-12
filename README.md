# 多基站传播图仿真

## 目录
- [简介](#简介)
- [使用方法](#使用方法)
- [依赖项](#依赖项)

## 简介
这个 MATLAB 文件实现了一个**多基站无线传播仿真与位置优化系统**，主要用途包括：

1. **加载校园/区域场景数据**（建筑物轮廓、初始基站位置等）。
2. **生成接收测试点**（避开建筑物内部）和**散射体**（墙体及随机分布，用于模拟多径传播）。
3. **通过粒子群算法（PSO）优化基站位置**，以最小化平均路径损耗、提高覆盖率和负载均衡度。
4. **计算传播损耗**：同时考虑直视路径（LOS）和散射多径路径（经过墙体/随机散射体），并利用三角形相交检测判断墙壁遮挡。
5. **结果可视化**：展示基站位置优化对比、PSO 收敛曲线、信号损耗分布图，并生成全区域覆盖热力图。
6. **保存所有结果**（优化位置、损耗矩阵、图片等）到指定目录。

该代码适用于**无线网络规划、室内外覆盖分析、基站选址优化**等科研或工程场景。

## 使用方法
你可以看到两个文件: **Multi_BS_Propagation_Graph_Simulation_Improved.m**和**summon.m**

**Multi_BS_Propagation_Graph_Simulation_Improved.m** 是主要的执行文件, 它读取指定路径下的基站&建筑文件信息并进行基站位置优化并生成示意图

**summon.m** 是地图资源文件生成器, 它在指定目录下生成一个*.mat*资源文件, 运行后可直接运行主执行文件来查看效果

## 依赖项
**Multi_BS_Propagation_Graph_Simulation_Improved.m** 运行必须满足下面的条件: 
1. 在指定位置(默认位置为D:\Propagation_code\Scenario_data\)存在data_TJSP.mat文件(可以修改如下代码块来改变读取文件路径):
```matlab
ROOT_PATH = 'D:\';
data_file = fullfile(ROOT_PATH, 'Propagation_code', 'Scenario_data', 'data_TJSP.mat');
if ~exist(data_file, 'file')
    error('❌ 未找到数据文件！\n   期望路径：%s\n   请先运行数据生成脚本。', data_file);
end
```
2. data_TJSP.mat必须符合以下格式:
GMsettings.Location_Tx : 基站初始位置 [N_bs x 3] (x,y,z)
XY                     : 建筑物轮廓点 [N_points x 2] (x,y)
building_start_index   : 每个建筑物在 XY 中的起始索引 [N_buildings x 1]
building_end_index     : 每个建筑物在 XY 中的结束索引 [N_buildings x 1]
num_of_buildings       : 建筑物数量
main_walls             : 墙体三角形列表 [N_walls x 12] (4个角点，每个角点 x,y,z)
可以直接修改summon.m中的参数来生成可用的文件:
### 1. 定义地图范围
```matlab
%% 1. 定义校园范围（示例：假设校园为 1000m x 800m 的矩形）
x_min = 0;   x_max = 1000;
y_min = 0;   y_max = 800;
```
### 2. 定义基站数目和位置
```matlab
%% 2. 生成模拟的基站初始位置（例如 3 个基站）
num_BS = 3;
GMsettings.Location_Tx = [
    200, 200, 10;    % 基站1 (x, y, z)
    500, 400, 10;    % 基站2
    800, 600, 10     % 基站3
];
```
### 3. 定义建筑物位置和尺寸
```matlab
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
```
### 4. 地图资源文件保存路径设置(可以更改名称)
```matlab
%% 5. 保存到 data_TJSP.mat
save_path = 'D:\Propagation_code\Scenario_data\data_TJSP.mat';
% 确保目标文件夹存在
folder = fileparts(save_path);
if ~exist(folder, 'dir')
    mkdir(folder);
 
end
```
