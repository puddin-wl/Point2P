function summary_root = run_mraf_wgs_yscale_probe()
% run_mraf_wgs_yscale_probe Tiny y-geometry calibration on the good WGS exp020 case.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = {
    struct('variant_name', "wgs_exp020_yscale1000", 'target_y_scale', 1.000)
    struct('variant_name', "wgs_exp020_yscale1015", 'target_y_scale', 1.015)
    struct('variant_name', "wgs_exp020_yscale1025", 'target_y_scale', 1.025)
    struct('variant_name', "wgs_exp020_yscale1035", 'target_y_scale', 1.035)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.variant_name, 'preset_name', "wgs_exp020_yscale_probe", 'artifact_group', "wgs_yscale_probe", ...
        'target_edge_mode', "pivot50", 'pivot50_inner_power', 0.75, 'pivot50_outer_power', 1.25, ...
        'target_x_scale', 1.0, 'target_y_scale', c.target_y_scale, ...
        'free_threshold', 0.135, 'edge_low_threshold', 0.135, ...
        'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
        'phase_blend', 0.25, 'method', "WGS-MRAF", 'use_wgs', true, 'wgs_exponent', 0.20);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {k-1, string(c.variant_name), c.target_y_scale, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, m.rms_90, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        fullprop.mraf_metrics.size50_x_um, fullprop.mraf_metrics.size50_y_um, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

t = cell2table(rows, 'VariableNames', {'case_index','variant_name','target_y_scale','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um','shoulder_peak_x','shoulder_peak_y','core_rms','rms_90', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core','efficiency_inside_signal','efficiency_inside_guard', ...
    'fullprop_size50_x_um','fullprop_size50_y_um','fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
t = add_yscale_status(t);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'wgs_yscale_probe', char(stamp + "_wgs_yscale_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'wgs_yscale_probe_metrics.csv'));
write_yscale_summary(fullfile(summary_root, 'wgs_yscale_probe_summary.txt'), t);
plot_yscale_summary(fullfile(summary_root, 'wgs_yscale_probe_summary_plot.png'), t);
plot_yscale_profile_overlay(fullfile(summary_root, 'wgs_yscale_profile_overlay.png'), t);
fprintf('\nWGS y-scale probe summary written to:\n%s\n', summary_root);
end

function t = add_yscale_status(t)
b = t(1,:);
status = strings(height(t),1); reason = strings(height(t),1);
for i=1:height(t)
    hard = strings(0,1); warn = strings(0,1);
    if abs(t.size50_y_um(i)-120) > 2.0, warn(end+1)="size50_y_not_120pm2"; end %#ok<AGROW>
    if abs(t.size50_x_um(i)-b.size50_x_um) > 3.0, hard(end+1)="size50_x_shift"; end %#ok<AGROW>
    if t.shoulder_peak_x(i) > b.shoulder_peak_x + 0.012, hard(end+1)="shoulder_x"; end %#ok<AGROW>
    if t.shoulder_peak_y(i) > b.shoulder_peak_y + 0.012, hard(end+1)="shoulder_y"; end %#ok<AGROW>
    if t.core_rms(i) > b.core_rms * 1.15, hard(end+1)="core_rms"; end %#ok<AGROW>
    if t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core + 0.04 || t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core + 0.04
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
t.delta_size50_y_to_120_um=t.size50_y_um-120;
t.delta_size50_y_vs_base_um=t.size50_y_um-b.size50_y_um;
t.delta_transition_x_um=t.transition_13p5_90_x_um-b.transition_13p5_90_x_um;
t.delta_transition_y_um=t.transition_13p5_90_y_um-b.transition_13p5_90_y_um;
t.delta_shoulder_x=t.shoulder_peak_x-b.shoulder_peak_x;
t.delta_shoulder_y=t.shoulder_peak_y-b.shoulder_peak_y;
t.delta_core_rms=t.core_rms-b.core_rms;
end

function write_yscale_summary(path,t)
fid=fopen(path,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'WGS y-scale probe summary\n=========================\n\n');
fprintf(fid,'Purpose: calibrate y-direction 50%% size toward 120 um without changing WGS exponent or pivot50 strength.\n');
fprintf(fid,'Fixed target: pivot50 p075/q125, free_threshold=0.135, n_iter=10, mraf_factor=1.0, phase_blend=0.25, WGS exponent=0.20.\n');
fprintf(fid,'Only target_y_scale is varied; RD baseline field is not geometrically scaled.\n\n');
for i=1:height(t)
    fprintf(fid,'%s target_y_scale=%.4f  %s\n',t.variant_name(i),t.target_y_scale(i),t.pass_fail_reason(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um, y-to-120 delta %.3f um; size13.5 %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.delta_size50_y_to_120_um(i),t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um, delta %.3f x %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i),t.delta_transition_x_um(i),t.delta_transition_y_um(i));
    fprintf(fid,'  shoulder %.6g x %.6g, delta %.6g x %.6g\n',t.shoulder_peak_x(i),t.shoulder_peak_y(i),t.delta_shoulder_x(i),t.delta_shoulder_y(i));
    fprintf(fid,'  core_rms %.6g, rms_90 %.6g, side_lobe %.6g x %.6g, efficiency %.6g / %.6g\n\n',t.core_rms(i),t.rms_90(i),t.side_lobe_peak_x_rel_to_core(i),t.side_lobe_peak_y_rel_to_core(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
ok=t(t.status~="FAIL",:);
if isempty(ok)
    fprintf(fid,'Suggested candidate: none; all y-scale cases failed hard constraints.\n');
else
    score=abs(ok.size50_y_um-120)+0.4*abs(ok.size50_x_um-330)+30*max(ok.delta_shoulder_x,0)+30*max(ok.delta_shoulder_y,0)+10*max(ok.delta_core_rms,0);
    [~,idx]=min(score); fprintf(fid,'Suggested candidate: %s\n',ok.variant_name(idx));
end
end

function plot_yscale_summary(path,t)
fig=figure('Visible','off','Color','w','Position',[100 100 1500 950]); tiledlayout(3,2,'TileSpacing','compact'); x=t.target_y_scale;
nexttile; plot_pair(x,t.size50_x_um,t.size50_y_um,'Size50 (um)','size50'); yline(330,'r--'); yline(120,'b--');
nexttile; plot_pair(x,t.transition_13p5_90_x_um,t.transition_13p5_90_y_um,'TW13.5-90 (um)','transition');
nexttile; plot_pair(x,t.shoulder_peak_x,t.shoulder_peak_y,'shoulder peak','shoulder');
nexttile; plot(x,t.core_rms,'ko-',x,t.rms_90,'mo-','LineWidth',1.2); grid on; xlabel('target y scale'); ylabel('RMS'); title('core RMS / rms90'); legend('core','rms90');
nexttile; plot_pair(x,t.side_lobe_peak_x_rel_to_core,t.side_lobe_peak_y_rel_to_core,'side lobe rel core','side lobe');
nexttile; plot(x,t.efficiency_inside_signal,'ro-',x,t.efficiency_inside_guard,'bo-','LineWidth',1.2); grid on; xlabel('target y scale'); ylabel('efficiency'); title('efficiency'); legend('signal','guard');
sgtitle('WGS exp020 y-scale probe: target_y_scale only'); exportgraphics(fig,path,'Resolution',150); close(fig);
end
function plot_pair(x,yx,yy,ylab,tit)
plot(x,yx,'ro-',x,yy,'bo-','LineWidth',1.2); grid on; xlabel('target y scale'); ylabel(ylab); title(tit); legend('x','y','Location','best'); xline(1.0,'k--');
end
function plot_yscale_profile_overlay(path,t)
fig=figure('Visible','off','Color','w','Position',[80 80 1400 820]); tiledlayout(1,2,'TileSpacing','compact'); ax1=nexttile; ax2=nexttile; colors=lines(height(t));
for row=1:height(t)
    xt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','x_profile_rd_target_mraf.csv')); yt=readtable(fullfile(string(t.artifact_root(row)),'ideal_fft','y_profile_rd_target_mraf.csv'));
    if row==1, plot(ax1,xt.x_m*1e6,xt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax1,'on'); plot(ax2,yt.y_m*1e6,yt.rd_intensity_core_norm,'k--','LineWidth',1.3,'DisplayName','RD baseline'); hold(ax2,'on'); end
    lw=1.0; if row==1, lw=1.8; end
    name=sprintf('yscale %.3f',t.target_y_scale(row));
    plot(ax1,xt.x_m*1e6,xt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',name);
    plot(ax2,yt.y_m*1e6,yt.mraf_intensity_core_norm,'-','Color',colors(row,:),'LineWidth',lw,'DisplayName',name);
end
format_axis(ax1,'X center profiles: WGS y-scale probe','x (um)',[-380 380]); format_axis(ax2,'Y center profiles: WGS y-scale probe','y (um)',[-230 230]);
sgtitle('WGS y-scale overlay: inspect 50% y crossing and shoulder'); exportgraphics(fig,path,'Resolution',150); close(fig);
end
function format_axis(ax,tit,xlab,xlims)
grid(ax,'on'); xlim(ax,xlims); ylim(ax,[0 1.35]); yline(ax,0.90,'g:'); yline(ax,0.50,'r:'); yline(ax,0.135,'m:'); title(ax,tit); xlabel(ax,xlab); ylabel(ax,'normalized intensity'); legend(ax,'Location','northeastoutside');
end
