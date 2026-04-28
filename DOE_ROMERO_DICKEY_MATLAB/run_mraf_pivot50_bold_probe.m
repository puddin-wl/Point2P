% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_pivot50_bold_probe()
% run_mraf_pivot50_bold_probe Directly test p075/q125 against p082/q118 baseline.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = {
    struct('variant_name', "pivot50_p082_q118_baseline", 'inner', 0.82, 'outer', 1.18)
    struct('variant_name', "pivot50_p075_q125", 'inner', 0.75, 'outer', 1.25)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.variant_name, 'preset_name', "caseC_stable", 'artifact_group', "pivot50_bold_probe", ...
        'target_edge_mode', "pivot50", 'pivot50_inner_power', c.inner, 'pivot50_outer_power', c.outer, ...
        'free_threshold', 0.135, 'edge_low_threshold', 0.135, ...
        'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
        'phase_blend', 0.25, 'method', "GS-MRAF", 'use_wgs', false, 'wgs_exponent', 0);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {k-1, string(c.variant_name), c.inner, c.outer, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

t = cell2table(rows, 'VariableNames', {'case_index','variant_name','inner_power','outer_power','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um','shoulder_peak_x','shoulder_peak_y','core_rms', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core','efficiency_inside_signal','efficiency_inside_guard', ...
    'fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
t = add_bold_status(t);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'pivot50_bold_probe', char(stamp + "_pivot50_bold_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'pivot50_bold_metrics.csv'));
write_bold_summary(fullfile(summary_root, 'pivot50_bold_summary.txt'), t);
plot_bold_overlay(fullfile(summary_root, 'pivot50_bold_profile_overlay.png'), t);
fprintf('\nPivot50 bold probe summary written to:\n%s\n', summary_root);
end

function t = add_bold_status(t)
b = t(1,:);
status = strings(height(t),1); reason = strings(height(t),1);
for i=1:height(t)
    hard = strings(0,1);
    if abs(t.size50_x_um(i)-b.size50_x_um) > 2, hard(end+1)="size50_x"; end %#ok<AGROW>
    if abs(t.size50_y_um(i)-b.size50_y_um) > 2, hard(end+1)="size50_y"; end %#ok<AGROW>
    if t.shoulder_peak_x(i) > b.shoulder_peak_x + 0.008, hard(end+1)="shoulder_x"; end %#ok<AGROW>
    if t.shoulder_peak_y(i) > b.shoulder_peak_y + 0.008, hard(end+1)="shoulder_y"; end %#ok<AGROW>
    if t.core_rms(i) > b.core_rms * 1.10, hard(end+1)="core_rms"; end %#ok<AGROW>
    if t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core + 0.03 || t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core + 0.03, hard(end+1)="side_lobe_large"; end %#ok<AGROW>
    if isempty(hard), status(i)="PASS"; reason(i)="PASS"; else, status(i)="FAIL"; reason(i)="FAIL: "+strjoin(hard,", "); end
end
t.status=status; t.pass_fail_reason=reason;
t.delta_size50_x_um=t.size50_x_um-b.size50_x_um; t.delta_size50_y_um=t.size50_y_um-b.size50_y_um;
t.delta_transition_x_um=t.transition_13p5_90_x_um-b.transition_13p5_90_x_um; t.delta_transition_y_um=t.transition_13p5_90_y_um-b.transition_13p5_90_y_um;
t.delta_shoulder_x=t.shoulder_peak_x-b.shoulder_peak_x; t.delta_shoulder_y=t.shoulder_peak_y-b.shoulder_peak_y;
t.delta_core_rms=t.core_rms-b.core_rms;
end

function write_bold_summary(path,t)
fid=fopen(path,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'Pivot50 bold p075/q125 probe\n==============================\n\n');
fprintf(fid,'Baseline row is p082/q118. Test row is p075/q125. Other params fixed: phase_blend=0.25, free_threshold=0.135, n_iter=10, gamma=1.0, mraf_factor=1.0, GS-MRAF, WGS off.\n\n');
for i=1:height(t)
    fprintf(fid,'%s inner=%.3f outer=%.3f  %s\n',t.variant_name(i),t.inner_power(i),t.outer_power(i),t.pass_fail_reason(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um, delta %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.delta_size50_x_um(i),t.delta_size50_y_um(i));
    fprintf(fid,'  size13.5 %.3f x %.3f um\n',t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um, delta %.3f x %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i),t.delta_transition_x_um(i),t.delta_transition_y_um(i));
    fprintf(fid,'  shoulder %.6g x %.6g, delta %.6g x %.6g\n',t.shoulder_peak_x(i),t.shoulder_peak_y(i),t.delta_shoulder_x(i),t.delta_shoulder_y(i));
    fprintf(fid,'  core_rms %.6g, side_lobe %.6g x %.6g, efficiency %.6g / %.6g\n\n',t.core_rms(i),t.side_lobe_peak_x_rel_to_core(i),t.side_lobe_peak_y_rel_to_core(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
end

function plot_bold_overlay(path,t)
fig=figure('Visible','off','Color','w','Position',[80 80 1350 800]); tiledlayout(1,2,'TileSpacing','compact'); ax1=nexttile; ax2=nexttile; colors=lines(height(t));
for row=1:height(t)
    xt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','x_profile_rd_target_mraf.csv')); yt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','y_profile_rd_target_mraf.csv'));
    if row==1, plot(ax1,xt.x_m*1e6,xt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax1,'on'); plot(ax2,yt.y_m*1e6,yt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax2,'on'); end
    lw=1.8; if row>1, lw=1.3; end
    plot(ax1,xt.x_m*1e6,xt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
    plot(ax2,yt.y_m*1e6,yt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
end
format_axis(ax1,'X center profiles: p082/q118 vs p075/q125','x (um)',[-380 380]); format_axis(ax2,'Y center profiles: p082/q118 vs p075/q125','y (um)',[-230 230]);
sgtitle('Pivot50 bold probe: inspect shoulder and side lobe'); exportgraphics(fig,path,'Resolution',150); close(fig);
end
function format_axis(ax,tit,xlab,xlims)
grid(ax,'on'); xlim(ax,xlims); ylim(ax,[0 1.35]); yline(ax,0.90,'g:'); yline(ax,0.50,'r:'); yline(ax,0.135,'m:'); title(ax,tit); xlabel(ax,xlab); ylabel(ax,'normalized intensity'); legend(ax,'Location','northeastoutside');
end
