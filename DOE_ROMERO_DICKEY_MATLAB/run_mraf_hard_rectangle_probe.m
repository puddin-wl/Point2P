function summary_root = run_mraf_hard_rectangle_probe()
% run_mraf_hard_rectangle_probe Hard-rectangle signal + MRAF free-region control experiment.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = {
    struct('variant_name', "baseline_rd_free0135", 'target_mode', "rd_derived", 'n_iter', 10)
    struct('variant_name', "hardrect_iter3", 'target_mode', "hard_rectangle", 'n_iter', 3)
    struct('variant_name', "hardrect_iter5", 'target_mode', "hard_rectangle", 'n_iter', 5)
    struct('variant_name', "hardrect_iter8", 'target_mode', "hard_rectangle", 'n_iter', 8)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.variant_name, 'preset_name', "caseC_stable", 'artifact_group', "hard_rectangle_probe", ...
        'target_mode', c.target_mode, 'target_edge_mode', "identity", ...
        'free_threshold', 0.135, 'edge_low_threshold', 0.135, ...
        'n_iter', c.n_iter, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
        'phase_blend', 0.25, 'method', "GS-MRAF", 'use_wgs', false, 'wgs_exponent', 0);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {k-1, string(c.variant_name), string(c.target_mode), c.n_iter, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, m.rms_90, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

t = cell2table(rows, 'VariableNames', {'case_index','variant_name','target_mode','n_iter','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um', ...
    'shoulder_peak_x','shoulder_peak_y','core_rms','rms_90', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core', ...
    'efficiency_inside_signal','efficiency_inside_guard', ...
    'fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
t = add_hardrect_deltas(t);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'hard_rectangle_probe', char(stamp + "_hardrect_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'hardrect_probe_metrics.csv'));
write_hardrect_summary(fullfile(summary_root, 'hardrect_probe_summary.txt'), t);
plot_hardrect_summary(fullfile(summary_root, 'hardrect_probe_summary_plot.png'), t);
plot_hardrect_profile_overlay(fullfile(summary_root, 'hardrect_profile_overlay.png'), t);
fprintf('\nHard rectangle probe summary written to:\n%s\n', summary_root);
end

function t = add_hardrect_deltas(t)
b = t(1,:);
t.delta_transition_x_um = t.transition_13p5_90_x_um - b.transition_13p5_90_x_um;
t.delta_transition_y_um = t.transition_13p5_90_y_um - b.transition_13p5_90_y_um;
t.delta_shoulder_x = t.shoulder_peak_x - b.shoulder_peak_x;
t.delta_shoulder_y = t.shoulder_peak_y - b.shoulder_peak_y;
t.delta_size50_x_um = t.size50_x_um - b.size50_x_um;
t.delta_size50_y_um = t.size50_y_um - b.size50_y_um;
t.delta_core_rms = t.core_rms - b.core_rms;
t.delta_side_lobe_x = t.side_lobe_peak_x_rel_to_core - b.side_lobe_peak_x_rel_to_core;
t.delta_side_lobe_y = t.side_lobe_peak_y_rel_to_core - b.side_lobe_peak_y_rel_to_core;
end

function write_hardrect_summary(path,t)
fid=fopen(path,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'Hard rectangle + few-step MRAF probe\n=====================================\n\n');
fprintf(fid,'Baseline: baseline_rd_free0135, rd_derived, free_threshold=0.135, n_iter=10.\n');
fprintf(fid,'Hard rectangle: signal |x|<=165 um, |y|<=60 um, amplitude=1; outside is MRAF free region (NaN), not zero-forced.\n');
fprintf(fid,'Fixed hardrect params: mraf_factor=1.0, phase_blend=0.25, GS-MRAF, WGS off.\n\n');
for i=1:height(t)
    fprintf(fid,'%s target=%s n_iter=%d\n',t.variant_name(i),t.target_mode(i),t.n_iter(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um, delta %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.delta_size50_x_um(i),t.delta_size50_y_um(i));
    fprintf(fid,'  size13.5 %.3f x %.3f um\n',t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um, delta %.3f x %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i),t.delta_transition_x_um(i),t.delta_transition_y_um(i));
    fprintf(fid,'  shoulder %.6g x %.6g, delta %.6g x %.6g\n',t.shoulder_peak_x(i),t.shoulder_peak_y(i),t.delta_shoulder_x(i),t.delta_shoulder_y(i));
    fprintf(fid,'  core_rms %.6g, rms_90 %.6g, side_lobe %.6g x %.6g, efficiency %.6g / %.6g\n\n',t.core_rms(i),t.rms_90(i),t.side_lobe_peak_x_rel_to_core(i),t.side_lobe_peak_y_rel_to_core(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
hard = t(t.target_mode=="hard_rectangle",:);
[~, idx] = min(hard.transition_13p5_90_x_um + hard.transition_13p5_90_y_um + 20*max(hard.delta_shoulder_x,0) + 20*max(hard.delta_shoulder_y,0) + 2*max(hard.delta_core_rms,0));
fprintf(fid,'Most reasonable hard-rectangle case by balanced score: %s\n',hard.variant_name(idx));
fprintf(fid,'Interpretation should consider profile overlay and ROI images; transition alone is not sufficient.\n');
end

function plot_hardrect_summary(path,t)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 1050]); tiledlayout(3,2,'TileSpacing','compact');
hard = t(t.target_mode=="hard_rectangle",:); x = hard.n_iter; b = t(1,:);
nexttile; plot_with_baseline(x, hard.transition_13p5_90_x_um, hard.transition_13p5_90_y_um, b.transition_13p5_90_x_um, b.transition_13p5_90_y_um, 'TW13.5-90 (um)', 'transition');
nexttile; plot_with_baseline(x, hard.shoulder_peak_x, hard.shoulder_peak_y, b.shoulder_peak_x, b.shoulder_peak_y, 'shoulder peak', 'shoulder');
nexttile; plot(x, hard.core_rms, 'ko-', 'LineWidth', 1.2); yline(b.core_rms,'k--','baseline'); grid on; xlabel('n iter'); ylabel('core RMS'); title('core RMS');
nexttile; plot_with_baseline(x, hard.size50_x_um, hard.size50_y_um, b.size50_x_um, b.size50_y_um, 'Size50 (um)', 'size50');
nexttile; plot_with_baseline(x, hard.side_lobe_peak_x_rel_to_core, hard.side_lobe_peak_y_rel_to_core, b.side_lobe_peak_x_rel_to_core, b.side_lobe_peak_y_rel_to_core, 'side lobe rel core', 'side lobe');
nexttile; plot(x, hard.efficiency_inside_signal, 'ro-', 'LineWidth', 1.2); yline(b.efficiency_inside_signal,'k--','baseline'); grid on; xlabel('n iter'); ylabel('efficiency inside signal'); title('efficiency');
sgtitle('Hard rectangle + few-step MRAF vs baseline_rd_free0135'); exportgraphics(fig,path,'Resolution',150); close(fig);
end

function plot_with_baseline(x,yx,yy,bx,by,ylab,tit)
plot(x,yx,'ro-',x,yy,'bo-','LineWidth',1.2); yline(bx,'r--','base x'); yline(by,'b--','base y'); grid on; xlabel('n iter'); ylabel(ylab); title(tit); legend('x','y','Location','best');
end

function plot_hardrect_profile_overlay(path,t)
fig=figure('Visible','off','Color','w','Position',[80 80 1350 800]); tiledlayout(1,2,'TileSpacing','compact'); ax1=nexttile; ax2=nexttile; colors=lines(height(t));
for row=1:height(t)
    xt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','x_profile_rd_target_mraf.csv')); yt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','y_profile_rd_target_mraf.csv'));
    if row==1
        plot(ax1,xt.x_m*1e6,xt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax1,'on');
        plot(ax2,yt.y_m*1e6,yt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax2,'on');
    end
    lw=1.0; if row==1, lw=1.8; end
    plot(ax1,xt.x_m*1e6,xt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
    plot(ax2,yt.y_m*1e6,yt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
end
format_axis(ax1,'X center profiles: RD, rd baseline, hard rectangle cases','x (um)',[-380 380]); format_axis(ax2,'Y center profiles: RD, rd baseline, hard rectangle cases','y (um)',[-230 230]);
sgtitle('Hard rectangle profile overlay: inspect transition, shoulder, side lobe'); exportgraphics(fig,path,'Resolution',150); close(fig);
end

function format_axis(ax,tit,xlab,xlims)
grid(ax,'on'); xlim(ax,xlims); ylim(ax,[0 1.6]); yline(ax,0.90,'g:'); yline(ax,0.50,'r:'); yline(ax,0.135,'m:'); title(ax,tit); xlabel(ax,xlab); ylabel(ax,'normalized intensity'); legend(ax,'Location','northeastoutside');
end