clear; close all; font_size=15;

file = 'temp.txt';

%根据MD模拟确定一些参数
num_blocks = 20;        % 块的个数
num_outputs = 10;       % 输出块温度的次数
Nx = 20;                % 输运方向晶胞个数
a = 5.4;                % 晶格常数，单位为 A
Lx = a * Nx;            % 输运方向长度
dx = Lx / num_blocks;   % 每个块的长度
power = 0.05 * 1.6e-7;  % 功率 eV/ps -> W
A = (4*a)^2 * 1.0e-20;  % 横截面积 A^2 -> m^2
Q = power / A / 2;      % 热流密度 W/m^2
t0 = 100;               % 每次输出数据增加的模拟时间，单位为  ps

% 非平衡模拟时间
t= t0 * (1 : num_outputs);

% 块坐标
x = (1 : num_blocks) * dx; % 一行
x = x.';                   % 一列

%从LAMMPS输出文件中提取块温度
[Chunk Coord1 Ncount v_TEMP] ...
    = textread(file, '%n%n%n%n', 'headerlines', 3, 'emptyvalue', 0);
temp = [Chunk Coord1 Ncount v_TEMP]; nall = length(temp) / (num_blocks+1);
clear Chunk Coord1 Ncount v_TEMP;
for i =1:nall
    temp((i-1)*num_blocks+1,:)=[];
end
temp = reshape(temp(:, end), num_blocks, num_outputs);

% 对输出次数作循环，计算温度梯度和热导率
for n = 1 : num_outputs
    
    % 确定拟合区间（具体问题具体分析）
    index_1 = 2 : num_blocks/2;              % 左边除去热源后的块指标
    index_2 = num_blocks/2+2 : num_blocks;   % 右边除去热汇后的块指标
    
    % 拟合温度分布
    p1 = fminsearch(@(p) ...
        norm(temp(index_1, n) - (p(1)-p(2)*x(index_1)) ), [60, -1]);
    p2 = fminsearch(@(p) ...
        norm(temp(index_2, n) - (p(1)+p(2)*x(index_2)) ), [60, 1]);
    
    % 得到温度梯度
    gradient_1 = p1(2) * 1.0e10; % K/A -> K/m
    gradient_2 = p2(2) * 1.0e10; % K/A -> K/m
    
    % 得到热导率
    kappa_1(n) = Q / gradient_1
    kappa_2(n) = Q / gradient_2
    
    % 画温度分布以及拟合曲线
    figure;
    plot(x, temp(:, n), 'bo', 'linewidth', 2);
    hold on;
    x1=x(index_1(1)) : 0.1 : x(index_1(end));
    x2=x(index_2(1)) : 0.1 : x(index_2(end));
    plot(x1, p1(1)-p1(2)*x1, 'r-', 'linewidth', 2);
    plot(x2, p2(1)+p2(2)*x2, 'g--', 'linewidth', 2);
    legend('MD', 'fit-left', 'fit-right');
    set(gca, 'fontsize', font_size);
    title (['t = ', num2str(t0 * n), ' ps']);
    xlabel('x (nm)');
    ylabel('T (K)');
end

% 热导率的时间收敛测试
kappa = (kappa_1 + kappa_2) / 2;
figure;
plot(t, kappa_1, ' rd', 'linewidth', 2);
hold on;
plot(t, kappa_2, ' bs', 'linewidth', 2);
plot(t, kappa, ' ko', 'linewidth', 2);
set(gca, 'fontsize', font_size);
xlabel('t (ps)');
ylabel('\kappa (W/mK)');
legend('left', 'right', 'average');

% 报道结果（具体问题具体分析）
kappa_average = mean( kappa(5:end) )
kappa_error   = mean( 0.5*abs( kappa_1(5:end)-kappa_2(5:end) ) )
