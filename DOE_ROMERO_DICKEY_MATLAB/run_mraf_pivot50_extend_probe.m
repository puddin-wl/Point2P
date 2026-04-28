function summary_root = run_mraf_pivot50_extend_probe()
% run_mraf_pivot50_extend_probe Extend pivot50 strength beyond p090/q110 conservatively.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = {
    struct('variant_name', "pivot50_p090_q110_baseline", 'inner', 0.90, 'outer', 1.10)
    struct('variant_name', "pivot50_p088_q112", 'inner', 0.88, 'outer', 1.12)
    struct('variant_name', "pivot50_p086_q114", 'inner', 0.86, 'outer', 1.14)
    struct('variant_name', "pivot50_p085_q115", 'inner', 0.85, 'outer', 1.15)
    struct('variant_name', "pivot50_p082_q118", 'inner', 0.82, 'outer', 1.18)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.variant_name, 'preset_name', "caseC_stable", 'artifact_group', "pivot50_extend_probe", ...
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
t = add_extend_status(t);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'pivot50_extend_probe', char(stamp + "_pivot50_extend_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'pivot50_extend_metrics.csv'));
write_extend_summary(fullfile(summary_root, 'pivot50_extend_summary.txt'), t);
plot_extend_summary(fullfile(summary_root, 'pivot50_extend_summary_plot.png'), t);
plot_extend_profile_overlay(fullfile(summary_root, 'pivot50_extend_profile_overlay.png'), t);
fprintf('\nPivot50 extend probe summary written to:\n%s\n', summary_root);
end

function t = add_extend_status(t)
b = t(1,:);
status = strings(height(t),1); reason = strings(height(t),1);
for i=1:height(t)
    hard = strings(0,1); warn = strings(0,1);
    if abs(t.size50_x_um(i)-b.size50_x_um) > 2, hard(end+1)="size50_x"; end %#ok<AGROW>
    if abs(t.size50_y_um(i)-b.size50_y_um) > 2, hard(end+1)="size50_y"; end %#ok<AGROW>
    if t.shoulder_peak_x(i) > b.shoulder_peak_x + 0.008, hard(end+1)="shoulder_x"; end %#ok<AGROW>
    if t.shoulder_peak_y(i) > b.shoulder_peak_y + 0.008, hard(end+1)="shoulder_y"; end %#ok<AGROW>
    if t.core_rms(i) > b.core_rms * 1.10, hard(end+1)="core_rms"; end %#ok<AGROW>
    if t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core + 0.03 || t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core + 0.03
        warn(end+1)="side_lobe"; %#ok<AGROW>
    end
    if t.efficiency_inside_signal(i) < b.efficiency_inside_signal * 0.98
        warn(end+1)="efficiency"; %#ok<AGROW>
    end
    if ~isempty(hard)
        status(i)="FAIL"; reason(i)="FAIL: "+strjoin(hard,", ");
    elseif ~isempty(warn)
        status(i)="PASS_WITH_WARNING"; reason(i)="WARN: "+strjoin(warn,", ");
    else
        status(i)="PASS"; reason(i)="PASS";
    end
end
t.status=status; t.pass_fail_reason=reason;
t.delta_size50_x_um=t.size50_x_um-b.size50_x_um; t.delta_size50_y_um=t.size50_y_um-b.size50_y_um;
t.delta_transition_x_um=t.transition_13p5_90_x_um-b.transition_13p5_90_x_um; t.delta_transition_y_um=t.transition_13p5_90_y_um-b.transition_13p5_90_y_um;
t.delta_shoulder_x=t.shoulder_peak_x-b.shoulder_peak_x; t.delta_shoulder_y=t.shoulder_peak_y-b.shoulder_peak_y;
t.delta_core_rms=t.core_rms-b.core_rms;
end

function write_extend_summary(path,t)
fid=fopen(path,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'Pivot50 strength extension summary\n==================================\n\n');
fprintf(fid,'Fixed: phase_blend=0.25, free_threshold=0.135, n_iter=10, gamma=1.0, mraf_factor=1.0, GS-MRAF, WGS off.\n');
fprintf(fid,'Baseline row is p090/q110. This intentionally tests stronger pivot50 remap only.\n\n');
for i=1:height(t)
    fprintf(fid,'%s inner=%.3f outer=%.3f  %s\n',t.variant_name(i),t.inner_power(i),t.outer_power(i),t.pass_fail_reason(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um, delta %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.delta_size50_x_um(i),t.delta_size50_y_um(i));
    fprintf(fid,'  size13.5 %.3f x %.3f um\n',t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um, delta %.3f x %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i),t.delta_transition_x_um(i),t.delta_transition_y_um(i));
    fprintf(fid,'  shoulder %.6g x %.6g, delta %.6g x %.6g\n',t.shoulder_peak_x(i),t.shoulder_peak_y(i),t.delta_shoulder_x(i),t.delta_shoulder_y(i));
    fprintf(fid,'  core_rms %.6g, side_lobe %.6g x %.6g, efficiency %.6g / %.6g\n\n',t.core_rms(i),t.side_lobe_peak_x_rel_to_core(i),t.side_lobe_peak_y_rel_to_core(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
ok=t(t.status~="FAIL",:);
if isempty(ok)
    fprintf(fid,'No stronger pivot50 candidate passed hard constraints; keep p090/q110.\n');
else
    score=ok.delta_transition_x_um+ok.delta_transition_y_um+60*max(ok.delta_shoulder_x,0)+60*max(ok.delta_shoulder_y,0)+0.3*abs(ok.delta_size50_x_um)+0.3*abs(ok.delta_size50_y_um);
    [~,idx]=min(score); fprintf(fid,'Suggested conservative candidate: %s\n',ok.variant_name(idx));
end
end

function plot_extend_summary(path,t)
fig=figure('Visible','off','Color','w','Position',[100 100 1400 1050]); tiledlayout(4,2,'TileSpacing','compact'); x=t.case_index;
nexttile; plot_pair(x,t.transition_13p5_90_x_um,t.transition_13p5_90_y_um,'TW13.5-90 (um)','transition');
nexttile; plot_pair(x,t.size50_x_um,t.size50_y_um,'Size50 (um)','size50');
nexttile; plot_pair(x,t.size13p5_x_um,t.size13p5_y_um,'Size13.5 (um)','size13.5');
nexttile; plot_pair(x,t.shoulder_peak_x,t.shoulder_peak_y,'shoulder peak','shoulder');
nexttile; plot(x,t.core_rms,'ko-','LineWidth',1.2); grid on; xlabel('case index'); ylabel('core RMS'); title('core RMS'); xline(0,'k--');
nexttile; plot_pair(x,t.side_lobe_peak_x_rel_to_core,t.side_lobe_peak_y_rel_to_core,'side lobe rel core','side lobe');
nexttile; plot(x,t.efficiency_inside_signal,'ro-',x,t.efficiency_inside_guard,'bo-','LineWidth',1.2); grid on; xlabel('case index'); ylabel('efficiency'); legend('signal','guard'); title('efficiency'); xline(0,'k--');
nexttile; axis off; text(0,1,status_text(t),'VerticalAlignment','top','FontName','Consolas');
sgtitle('Pivot50 strength extension beyond p090/q110'); exportgraphics(fig,path,'Resolution',150); close(fig);
end
function plot_pair(x,yx,yy,ylab,tit)
plot(x,yx,'ro-',x,yy,'bo-','LineWidth',1.2); grid on; xlabel('case index'); ylabel(ylab); title(tit); legend('x','y','Location','best'); xline(0,'k--');
end
function txt=status_text(t)
txt="idx variant status"+newline; for i=1:height(t), txt=txt+sprintf('%d %s %s\n',t.case_index(i),t.variant_name(i),t.status(i)); end
end
function plot_extend_profile_overlay(path,t)
fig=figure('Visible','off','Color','w','Position',[80 80 1350 800]); tiledlayout(1,2,'TileSpacing','compact'); ax1=nexttile; ax2=nexttile; colors=lines(height(t));
for row=1:height(t)
    xt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','x_profile_rd_target_mraf.csv')); yt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','y_profile_rd_target_mraf.csv'));
    if row==1, plot(ax1,xt.x_m*1e6,xt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax1,'on'); plot(ax2,yt.y_m*1e6,yt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax2,'on'); end
    lw=1.0; if row==1, lw=1.8; end
    plot(ax1,xt.x_m*1e6,xt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
    plot(ax2,yt.y_m*1e6,yt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',char(t.variant_name(row)));
end
format_axis(ax1,'X center profiles: pivot50 strength extension','x (um)',[-380 380]); format_axis(ax2,'Y center profiles: pivot50 strength extension','y (um)',[-230 230]);
sgtitle('Pivot50 extension overlay: inspect size50, 13.5%, shoulder'); exportgraphics(fig,path,'Resolution',150); close(fig);
end
function format_axis(ax,tit,xlab,xlims)
grid(ax,'on'); xlim(ax,xlims); ylim(ax,[0 1.35]); yline(ax,0.90,'g:'); yline(ax,0.50,'r:'); yline(ax,0.135,'m:'); title(ax,tit); xlabel(ax,xlab); ylabel(ax,'normalized intensity'); legend(ax,'Location','northeastoutside');
end