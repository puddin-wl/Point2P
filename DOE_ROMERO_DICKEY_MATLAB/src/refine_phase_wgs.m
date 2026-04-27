function phase_refined_rad = refine_phase_wgs(phase_initial_rad, varargin)
% refine_phase_wgs 预留的轻量 WGS 相位 refinement 入口。
%
% 第一版项目默认不启用优化器。后续若加入少量迭代，应以解析相位为初值，
% 并始终保存原始 Romero-Dickey 结果用于对比。

phase_refined_rad = phase_initial_rad;
warning('refine_phase_wgs:NotImplemented', 'WGS refinement is reserved but not implemented/enabled in this first version.');
end

