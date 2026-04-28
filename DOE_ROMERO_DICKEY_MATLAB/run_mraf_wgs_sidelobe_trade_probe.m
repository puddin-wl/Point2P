% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_wgs_sidelobe_trade_probe()
% run_mraf_wgs_sidelobe_trade_probe Allow limited side-lobe budget to reduce shoulder_y.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = {
    struct('variant_name', "sidelobe_trade_baseline", 'y_free', false, 'half_y', 72, 'max_I', 0.50, 'damp', false, 'damp_factor', 1.00)
    struct('variant_name', "sidelobe_trade_A1_yfree_mild", 'y_free', true, 'half_y', 72, 'max_I', 0.50, 'damp', false, 'damp_factor', 1.00)
    struct('variant_name', "sidelobe_trade_A2_yfree_stronger", 'y_free', true, 'half_y', 68, 'max_I', 0.60, 'damp', false, 'damp_factor', 1.00)
    struct('variant_name', "sidelobe_trade_B1_yedge_damp", 'y_free', false, 'half_y', 72, 'max_I', 0.50, 'damp', true, 'damp_factor', 0.50)
    struct('variant_name', "sidelobe_trade_B2_yfree_damp", 'y_free', true, 'half_y', 72, 'max_I', 0.50, 'damp', true, 'damp_factor', 0.50)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.variant_name, 'preset_name', "wgs_sidelobe_trade", 'artifact_group', "wgs_sidelobe_trade_probe", ...
        'target_edge_mode', "pivot50", 'pivot50_inner_power', 0.75, 'pivot50_outer_power', 1.25, ...
        'target_x_scale', 1.0, 'target_y_scale', 1.022, ...
        'free_threshold', 0.135, 'edge_low_threshold', 0.135, ...
        'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
        'phase_blend', 0.25, 'method', "WGS-MRAF", 'use_wgs', true, 'wgs_exponent', 0.20, ...
        'y_free_reservoir_enable', c.y_free, 'y_free_reservoir_half_y_um', c.half_y, 'y_free_reservoir_max_intensity', c.max_I, ...
        'wgs_y_edge_damping_enable', c.damp, 'wgs_y_edge_damping_factor', c.damp_factor, ...
        'wgs_y_edge_damping_low', 0.50, 'wgs_y_edge_damping_high', 0.90);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {k-1, string(c.variant_name), c.y_free, c.half_y, c.max_I, c.damp, c.damp_factor, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, m.rms_90, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        fullprop.mraf_metrics.size50_x_um, fullprop.mraf_metrics.size50_y_um, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

t = cell2table(rows, 'VariableNames', {'case_index','variant_name','y_free_reservoir','y_free_half_y_um','y_free_max_I','wgs_y_edge_damping','wgs_y_edge_damping_factor','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um','shoulder_peak_x','shoulder_peak_y','core_rms','rms_90', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core','efficiency_inside_signal','efficiency_inside_guard', ...
    'fullprop_size50_x_um','fullprop_size50_y_um','fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
t = add_sidelobe_trade_status(t);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'wgs_sidelobe_trade_probe', char(stamp + "_sidelobe_trade_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'sidelobe_trade_metrics.csv'));
write_sidelobe_trade_summary(fullfile(summary_root, 'sidelobe_trade_summary.txt'), t);
plot_sidelobe_trade_summary(fullfile(summary_root, 'sidelobe_trade_summary_plot.png'), t);
plot_sidelobe_trade_profile_overlay(fullfile(summary_root, 'sidelobe_trade_profile_overlay.png'), t);
fprintf('\nSide-lobe trade probe summary written to:\n%s\n', summary_root);
end

function t = add_sidelobe_trade_status(t)
b = t(1,:);
status = strings(height(t),1); reason = strings(height(t),1);
for i=1:height(t)
    hard = strings(0,1); warn = strings(0,1);
    if t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core + 0.03, hard(end+1)="side_lobe_y_budget"; end %#ok<AGROW>
    if abs(t.size50_y_um(i)-120) > 2.0, hard(end+1)="size50_y"; end %#ok<AGROW>
    if abs(t.size50_x_um(i)-b.size50_x_um) > 2.0, hard(end+1)="size50_x_shift"; end %#ok<AGROW>
    if t.core_rms(i) > b.core_rms * 1.10, hard(end+1)="core_rms"; end %#ok<AGROW>
    if t.shoulder_peak_y(i) > b.shoulder_peak_y, warn(end+1)="shoulder_y_not_lower"; end %#ok<AGROW>
    if t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core + 0.03, warn(end+1)="side_lobe_x"; end %#ok<AGROW>
    if t.efficiency_inside_signal(i) < b.efficiency_inside_signal * 0.98, warn(end+1)="efficiency"; end %#ok<AGROW>
    if ~isempty(hard)
        status(i)="FAIL"; reason(i)="FAIL: "+strjoin(hard,", ");
    elseif ~isempty(warn)
        status(i)="PASS_WITH_WARNING"; reason(i)="WARN: "+strjoin(warn,", ");
    else
        status(i)="PASS"; reason(i)="PASS";
    end
end
t.status=status; t.pass_fail_reason=reason;
t.delta_shoulder_y=t.shoulder_peak_y-b.shoulder_peak_y;
t.delta_side_lobe_y=t.side_lobe_peak_y_rel_to_core-b.side_lobe_peak_y_rel_to_core;
t.delta_size50_y_um=t.size50_y_um-b.size50_y_um;
t.delta_transition_y_um=t.transition_13p5_90_y_um-b.transition_13p5_90_y_um;
t.delta_core_rms=t.core_rms-b.core_rms;
end

function write_sidelobe_trade_summary(path,t)
fid=fopen(path,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'WGS side-lobe trade probe summary\n==================================\n\n');
fprintf(fid,'Baseline: pivot50 p075/q125, free_threshold=0.135, target_y_scale=1.022, WGS exponent=0.20.\n');
fprintf(fid,'Goal: allow side_lobe_y up to baseline +0.03 to lower shoulder_y.\n\n');
for i=1:height(t)
    fprintf(fid,'%s  %s\n',t.variant_name(i),t.pass_fail_reason(i));
    fprintf(fid,'  knobs: y_free=%d half_y=%.1f um max_I=%.2f; damping=%d factor=%.2f\n',t.y_free_reservoir(i),t.y_free_half_y_um(i),t.y_free_max_I(i),t.wgs_y_edge_damping(i),t.wgs_y_edge_damping_factor(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um; size13.5 %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um, delta_y %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i),t.delta_transition_y_um(i));
    fprintf(fid,'  shoulder_y %.6g, delta %.6g; side_lobe_y %.6g, delta %.6g\n',t.shoulder_peak_y(i),t.delta_shoulder_y(i),t.side_lobe_peak_y_rel_to_core(i),t.delta_side_lobe_y(i));
    fprintf(fid,'  shoulder_x %.6g; side_lobe_x %.6g; core_rms %.6g; efficiency %.6g / %.6g\n\n',t.shoulder_peak_x(i),t.side_lobe_peak_x_rel_to_core(i),t.core_rms(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
ok=t(t.status~="FAIL",:);
if isempty(ok)
    fprintf(fid,'Suggested candidate: none; all cases failed hard constraints.\n');
else
    score=ok.delta_shoulder_y + 0.4*max(ok.delta_side_lobe_y,0) + 8*max(ok.delta_core_rms,0) + 0.05*abs(ok.delta_size50_y_um);
    [~,idx]=min(score);
    fprintf(fid,'Suggested candidate: %s\n',ok.variant_name(idx));
end
end

function plot_sidelobe_trade_summary(path,t)
fig=figure('Visible','off','Color','w','Position',[100 100 1500 950]); tiledlayout(3,2,'TileSpacing','compact'); x=t.case_index;
labels = cellstr(t.variant_name);
nexttile; plot(x,t.shoulder_peak_y,'bo-',x,t.shoulder_peak_x,'ro-','LineWidth',1.2); grid on; title('shoulder peak'); ylabel('shoulder'); legend('y','x'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
nexttile; plot(x,t.side_lobe_peak_y_rel_to_core,'bo-',x,t.side_lobe_peak_x_rel_to_core,'ro-','LineWidth',1.2); grid on; title('side lobe rel core'); ylabel('side lobe'); legend('y','x'); yline(t.side_lobe_peak_y_rel_to_core(1)+0.03,'b--','y budget'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
nexttile; plot(x,t.size50_y_um,'bo-',x,t.size50_x_um,'ro-','LineWidth',1.2); grid on; title('size50'); ylabel('um'); yline(120,'b--'); yline(330,'r--'); legend('y','x'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
nexttile; plot(x,t.core_rms,'ko-',x,t.rms_90,'mo-','LineWidth',1.2); grid on; title('RMS'); ylabel('RMS'); legend('core','rms90'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
nexttile; plot(x,t.transition_13p5_90_y_um,'bo-',x,t.transition_13p5_90_x_um,'ro-','LineWidth',1.2); grid on; title('TW13.5-90'); ylabel('um'); legend('y','x'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
nexttile; plot(x,t.efficiency_inside_signal,'ro-',x,t.efficiency_inside_guard,'bo-','LineWidth',1.2); grid on; title('efficiency'); ylabel('efficiency'); legend('signal','guard'); set(gca,'XTick',x,'XTickLabel',labels); xtickangle(25);
sgtitle('WGS side-lobe trade probe'); exportgraphics(fig,path,'Resolution',150); close(fig);
end

function plot_sidelobe_trade_profile_overlay(path,t)
fig=figure('Visible','off','Color','w','Position',[80 80 1400 820]); tiledlayout(1,2,'TileSpacing','compact'); ax1=nexttile; ax2=nexttile; colors=lines(height(t));
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
format_axis(ax1,'X center profiles: side-lobe trade','x (um)',[-380 380]);
format_axis(ax2,'Y center profiles: side-lobe trade','y (um)',[-230 230]);
sgtitle('Side-lobe trade overlay: inspect shoulder_y and outer lobes'); exportgraphics(fig,path,'Resolution',150); close(fig);
end

function format_axis(ax,tit,xlab,xlims)
grid(ax,'on'); xlim(ax,xlims); ylim(ax,[0 1.35]); yline(ax,0.90,'g:'); yline(ax,0.50,'r:'); yline(ax,0.135,'m:'); title(ax,tit); xlabel(ax,xlab); ylabel(ax,'normalized intensity'); legend(ax,'Location','northeastoutside');
end

