
% 清空工作区，关闭所有图形
clear; close all; clc;

% 提示用户输入文件名
filename = input('outputfile317_jizhibaohe：', 's');

% 检查文件是否存在
if ~exist(filename, 'file')
    error('文件 %s 不存在，请检查路径和文件名。', filename);
end

% 打开文件
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件 %s。', filename);
end

% 初始化存储数值的数组
valuesX = [];
valuesY = [];
values_gg = [];

% 逐行读取文件
while ~feof(fid)
    line = fgetl(fid);  % 读取一行
    if ischar(line) && contains(line, 'Current in x direction =')  % 检查是否包含目标字符串
        % 使用正则表达式提取等号后的数值（支持整数、小数、科学计数法）
        tokens = regexp(line, 'Current in x direction = \s*([-+]?\d*\.?\d+([eE][-+]?\d+)?)', 'tokens');
        if ~isempty(tokens)
            numStr = tokens{1}{1};  % 提取匹配的数字字符串
            num = str2double(numStr);  % 转换为数值
            if ~isnan(num)  % 确保转换成功
                valuesX = [valuesX; num];  % 添加到数组
            end
        end
    end
        % 提取 "gg =" 后的数值
    if ischar(line) && contains(line, 'Current in y direction =')  % 检查是否包含目标字符串
        % 使用正则表达式提取等号后的数值（支持整数、小数、科学计数法）
        tokens = regexp(line, 'Current in y direction = \s*([-+]?\d*\.?\d+([eE][-+]?\d+)?)', 'tokens');
        if ~isempty(tokens)
            numStr = tokens{1}{1};  % 提取匹配的数字字符串
            num = str2double(numStr);  % 转换为数值
            if ~isnan(num)  % 确保转换成功
                valuesY = [valuesY; num];  % 添加到数组
            end
        end
    end
                    % 提取 "gg =" 后的数值
     if ischar(line) && contains(line, 'gg =')  % 检查是否包含目标字符串
        if contains(line, 'gg =')
            tokens = regexp(line, 'gg =\s*([-+]?\d*\.?\d+([eE][-+]?\d+)?)', 'tokens');
            if ~isempty(tokens)
                numStr = tokens{1}{1};
                num = str2double(numStr);
                if ~isnan(num)
                    values_gg = [values_gg; num];
                end
            end
        end
     end
end

% 关闭文件
fclose(fid);

% 检查是否提取到任何数值
if isempty(valuesX)
    error('文件中没有找到任何 "Current in x direction " 的数值。');
end

% 将数值保存到文本文件，每行一个数字
outputFilename = 'extracted_values.txt';
dlmwrite(outputFilename, valuesX, 'precision', '%.6f');  % 保留6位小数，可根据需要调整
fprintf('已提取 %d 个数值，并保存到 %s\n', length(valuesX), outputFilename);

% 绘制数值变化图
figure;
plot(valuesX, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('数据点序号');
ylabel('数值');
title('x current 数值变化规律');
grid on;

figure;
plot(valuesY, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('数据点序号');
ylabel('数值');
title('y current 数值变化规律');
grid on;

figure();
semilogy(values_gg, 'r', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('数据点序号');
ylabel('数值');
title('gg 数值变化规律');
grid on


% 可选：如果需要绘制连续折线，也可使用 stem 图等
% stem(values, 'filled'); 

disp('绘图完成。');