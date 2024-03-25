function [ ] = Expand_axis_fill_figure( axis_handle )  %函数定义
% TightInset的位置
inset_vectior = get(axis_handle, 'TightInset');
% inset_x = 0.0174;%inset_vectior(1);
% inset_y = 0.0153;%inset_vectior(2);
% inset_w = 0.0046;%inset_vectior(3);
% inset_h = 0.0032;%inset_vectior(4);
inset_x = inset_vectior(1);
inset_y = inset_vectior(2);
inset_w = inset_vectior(3);
inset_h = inset_vectior(4);

% inset_x = inset_vectior(1)-0.03;
% inset_y = inset_vectior(2)-0.02;
% inset_w = inset_vectior(3)-0.005;
% inset_h = inset_vectior(4)-0.005;


% OuterPosition的位置
outer_vector = get(axis_handle, 'OuterPosition');
pos_new_x = outer_vector(1) + inset_x; % 将Position的原点移到到TightInset的原点
pos_new_y = outer_vector(2) + inset_y;
pos_new_w = outer_vector(3) - inset_w - inset_x; % 重设Position的宽
pos_new_h = outer_vector(4) - inset_h - inset_y; % 重设Position的高

% 重设Position
set(axis_handle, 'Position', [pos_new_x, pos_new_y, pos_new_w, pos_new_h]);
%函数结束