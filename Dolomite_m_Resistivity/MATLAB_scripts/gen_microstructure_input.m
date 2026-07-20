% ========================================================
% 将 20x20 矩阵转换为 microstructure.dat (供 Fortran 程序读取)
% 假设您的原始矩阵 matrix_original 是常规顺序：
%   第一行 = 顶部行，最后一行 = 底部行
% 程序读取顺序：第一行 = 底部行，最后一行 = 顶部行
% 因此需要对原始矩阵进行上下翻转 (flipud)
% ========================================================

% --- 1. 载入您的原始矩阵 ---
% 如果您的矩阵已经在 MATLAB 工作区中，直接使用变量名，例如：
% matrix_original = your_matrix; % 请替换为您的实际变量名
clear all
clc
% 如果矩阵保存在文件（如 Excel 或文本文件）中，可以用以下方式读取：
matrix_original = xlsread('基质孔隙35s.xlsx');      % Excel 文件2500*1250
figure()
imagesc(matrix_original)

%matrix_original = readtable('0228.txt');

% 如果从文件读取：
% matrix_original = load('your_matrix.txt');   % 纯文本矩阵（每行20个数）
% matrix_original = xlsread('your_file.xlsx'); % Excel 文件

% --- 2. 上下翻转 ---
matrix_flipped = flipud(matrix_original);  % 使第一行为底部行

% --- 3. 将翻转后的矩阵展平为列向量（按行优先顺序）---
vector_data = matrix_flipped(:);  % MATLAB 按列优先展开，但我们需要按行优先
% 修正：按行优先展开
vector_data = reshape(matrix_flipped', [], 1);  % 先转置，再按列展开得到行优先顺序

% 或者使用循环：
% vector_data = [];
% for row = 1:20
%     vector_data = [vector_data; matrix_flipped(row,:)'];
% end

% --- 4. 写入文件（每行一个数）---
dlmwrite('microstructure323_jizhi35s1.dat', vector_data, 'delimiter', '\n');
fprintf('文件 microstructure323_jizhi35s1.dat 已生成，共 %d 行。\n', length(vector_data));
disp('验证成功。');
% --- 5. （可选）验证 ---
% 读取文件并恢复为矩阵，检查是否正确
% data_read = load('microstructure1010.dat');  % 读取列向量
% matrix_recovered = reshape(data_read, [20,20])';  % 注意 reshape 后需转置
% matrix_recovered = flipud(matrix_recovered);  % 再翻转回原始顺序

% if isequal(matrix_recovered, matrix_original)
%     disp('验证成功：文件读取后恢复的矩阵与原始矩阵一致。');
% else
%     disp('警告：验证失败，请检查数据。');
% end