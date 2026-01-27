function set_XY_Tick

set(gca, 'Box', 'off', ...                                        % 边框
         'LineWidth', 1,...                                       % 线宽
         'XGrid', 'off', 'YGrid', 'off', ...                      % 网格
         'TickDir', 'out', 'TickLength', [.015 .015], ...         % 刻度
         'XMinorTick', 'on', 'YMinorTick', 'on', ...              % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])             % 坐标轴颜色