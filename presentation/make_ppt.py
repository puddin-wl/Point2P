"""Generate project presentation PPT — physics-focused, equations as rendered images."""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from io import BytesIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTPUT = "E:/program/Point2P/presentation/Point2P_Project_Report.pptx"
FONT = "Microsoft YaHei"
BLUE = RGBColor(0x1A, 0x56, 0x8E)
DARK = RGBColor(0x2D, 0x2D, 0x2D)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_BG = RGBColor(0xF5, 0xF7, 0xFA)
ACCENT = RGBColor(0xE8, 0x6A, 0x17)
GRAY = RGBColor(0x66, 0x66, 0x66)

# Make matplotlib use a math font that renders well
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.size"] = 12


def render_equation(latex, fontsize=14, dpi=150, color="black", pad=0.15):
    """Render a LaTeX equation string to a PNG image in memory."""
    fig, ax = plt.subplots(figsize=(0.01, 0.01))
    text = ax.text(0.5, 0.5, f"${latex}$", fontsize=fontsize, ha="center", va="center",
                   transform=ax.transAxes, color=color)
    ax.axis("off")
    fig.canvas.draw()
    bbox = text.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    plt.close(fig)
    w, h = bbox.width, bbox.height
    fig, ax = plt.subplots(figsize=(w + 2 * pad, h + 2 * pad), dpi=dpi)
    ax.text(0.5, 0.5, f"${latex}$", fontsize=fontsize, ha="center", va="center",
            transform=ax.transAxes, color=color)
    ax.axis("off")
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                transparent=True)
    plt.close(fig)
    buf.seek(0)
    return buf


def add_eq_image(slide, latex, left, top, width=None, fontsize=14, dpi=150):
    """Add a rendered equation as an inline image."""
    buf = render_equation(latex, fontsize=fontsize, dpi=dpi)
    if width is None:
        slide.shapes.add_picture(buf, Inches(left), Inches(top))
    else:
        slide.shapes.add_picture(buf, Inches(left), Inches(top), Inches(width))


def slide_layout(prs, idx=5):
    return prs.slide_layouts[idx]


def add_text(slide, text, left=0.6, top=1.1, width=8.8, height=5.5, size=15, bold=False, color=DARK):
    txBox = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, line in enumerate(text.strip().split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line
        p.font.size = Pt(size)
        p.font.name = FONT
        p.font.color.rgb = color
        p.font.bold = bold
        p.space_after = Pt(5)
    return txBox


def add_table(slide, headers, rows, left=0.6, top=1.5, width=8.8, row_height=0.40,
              col_widths=None, font_size=12):
    n_rows = len(rows) + 1
    n_cols = len(headers)
    table_shape = slide.shapes.add_table(n_rows, n_cols, Inches(left), Inches(top),
                                         Inches(width), Inches(row_height * n_rows))
    table = table_shape.table
    if col_widths:
        for i, w in enumerate(col_widths):
            table.columns[i].width = Inches(w)
    for j, h in enumerate(headers):
        cell = table.cell(0, j)
        cell.text = h
        cell.fill.solid()
        cell.fill.fore_color.rgb = BLUE
        for p in cell.text_frame.paragraphs:
            p.font.size = Pt(font_size)
            p.font.bold = True
            p.font.color.rgb = WHITE
            p.font.name = FONT
            p.alignment = PP_ALIGN.CENTER
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = table.cell(i + 1, j)
            cell.text = str(val)
            if i % 2 == 0:
                cell.fill.solid()
                cell.fill.fore_color.rgb = LIGHT_BG
            for p in cell.text_frame.paragraphs:
                p.font.size = Pt(font_size)
                p.font.name = FONT
                p.font.color.rgb = DARK
                p.alignment = PP_ALIGN.CENTER
    return table_shape


def add_bullet(slide, items, left=0.6, top=1.15, width=8.8, size=14, spacing=6):
    """Add bullet-point text."""
    lines = []
    for item in items:
        lines.append("• " + item)
    add_text(slide, "\n".join(lines), left=left, top=top, width=width, size=size)
    # Adjust space after each line
    # (python-pptx basic; bullets are just text here)


def make_ppt():
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    # ====================================================================
    # Slide 1: Cover
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    bg = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width, prs.slide_height)
    bg.fill.solid(); bg.fill.fore_color.rgb = BLUE; bg.line.fill.background()
    add_text(slide, "矩形平顶 DOE 光束整形\n—— 从设计到容差评估", left=0.8, top=1.8, width=8.5, height=2,
             size=36, bold=True, color=WHITE)
    add_text(slide, "基于 Romero-Dickey 点对点法 + RTAD 目标 + WGS 优化\n\n2026 年 5 月", left=0.8, top=4.0, width=8.5, height=1.5,
             size=18, color=WHITE)

    # ====================================================================
    # Slide 2: Project Goal
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "项目目标", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, "设计一片衍射光学元件（DOE），将 532nm 高斯激光束整形成焦平面上的矩形平顶光斑。", top=0.85, size=14, bold=True, color=ACCENT)
    add_table(slide,
              ["参数", "符号", "值"],
              [
                  ["激光波长", "λ", "532 nm"],
                  ["入射光斑 (1/e² 强度直径)", "2w₀", "5 mm"],
                  ["傅里叶透镜焦距", "f", "429 mm"],
                  ["通光孔径 (直径)", "D", "15 mm"],
                  ["目标光斑 (50% 全宽)", "W₅₀ × H₅₀", "330 × 120 μm"],
                  ["计算网格", "N × N", "2048 × 2048"],
                  ["焦面采样", "dx_focal", "2.5 μm/pixel"],
                  ["Romero-Dickey βₓ / βⵧ", "β", "9.06 / 3.30"],
              ],
              col_widths=[3.0, 2.0, 3.5], font_size=12, top=1.8, row_height=0.42)
    add_text(slide, "核心流程：  Point2P 初始相位  →  WGS 迭代优化  →  真实光路容差评估", top=6.5, size=14, bold=True, color=ACCENT)

    # ====================================================================
    # Slide 3: Romero-Dickey — concept + β
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ① — Romero-Dickey 点对点法", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, "Romero & Dickey, \"Lossless laser beam shaping,\" J. Opt. Soc. Am. A, 13(4), 751–760 (1996)", top=0.85, size=10, color=GRAY)
    add_bullet(slide, [
        "菲涅尔近似下，用稳相法 (stationary phase) 解析求解高斯→平顶的 DOE 相位分布",
        "解的质量由无量纲参数 β 决定：",
    ], top=1.1, size=14)

    add_eq_image(slide, r"\beta = \frac{2\pi \cdot r_i \cdot R_o}{\lambda \cdot f}",
                 left=0.8, top=2.0, fontsize=16, dpi=180)

    add_bullet(slide, [
        "β 越大 → 几何光学近似越精确，平顶边缘越陡峭（裙边宽度 ∝ 1/β⁰·⁹）",
        "β < 6/π ≈ 1.91 时无合理的平顶解（不确定性关系限制）",
    ], top=2.7, size=14)

    add_table(slide, ["参数", "符号", "值", "说明"],
              [
                  ["输入 1/e 幅度半径", "rᵢ", "1.768 mm", "= (5mm/2) / √2"],
                  ["x 方向输出尺度", "Roₓ", "186.2 μm", "330 μm / √π"],
                  ["y 方向输出尺度", "Roⵧ", "67.7 μm", "120 μm / √π"],
                  ["x 方向 β", "βₓ", "9.06", "正常"],
                  ["y 方向 β", "βⵧ", "3.30", "限制方向 (βⵧ < βₓ)"],
              ],
              col_widths=[2.3, 1.0, 1.8, 3.7], font_size=11, top=4.0, row_height=0.38)

    # ====================================================================
    # Slide 4: Romero-Dickey — phase formula
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ① — Romero-Dickey 相位公式", top=0.3, size=28, bold=True, color=BLUE)
    add_bullet(slide, [
        "可分离假设: φ(x,y) = φₓ(x) + φⵧ(y)，x 和 y 方向独立求解",
    ], top=1.0, size=14)

    add_text(slide, "一维无纲量相位 (xi = x / rᵢ)：", top=1.6, size=14, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"\phi(\xi) = \beta \cdot \left[ \xi \cdot \frac{\sqrt{\pi}}{2} \cdot \mathrm{erf}(\xi) "
        r"+ \frac{1}{2} \cdot e^{-\xi^{2}} - \frac{1}{2} \right]",
        left=0.5, top=2.0, fontsize=15, dpi=180)

    add_text(slide, "其中 erf(ξ) 为误差函数:  erf(ξ) = (2/√π) ∫₀^ξ exp(-t²) dt", top=2.7, size=12, color=GRAY)
    add_text(slide, "二维可分离相位：", top=3.3, size=14, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"\phi(x,y) = \beta_x \cdot f\!\left(\frac{x}{r_i}\right) \;+\; \beta_y \cdot f\!\left(\frac{y}{r_i}\right),"
        r"\quad f(\xi) = \xi\cdot\frac{\sqrt{\pi}}{2}\,\mathrm{erf}(\xi) + \frac{1}{2}e^{-\xi^{2}} - \frac{1}{2}",
        left=0.3, top=3.7, fontsize=14, dpi=180)

    add_bullet(slide, [
        "该相位为解析解，无需迭代，直接给出接近目标能量分布的初始相位 (phase0)",
    ], top=4.8, size=14)

    # ====================================================================
    # Slide 5: Romero-Dickey — sensitivity
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ① — 离焦与偏移敏感度 (Romero-Dickey)", top=0.3, size=26, bold=True, color=BLUE)
    add_bullet(slide, [
        "Romero & Dickey 在论文中也推导了系统误差的敏感度分析，本项目在真实光路模拟中独立验证",
    ], top=0.95, size=13)

    add_text(slide, "轴向离焦 (dz)：", top=1.5, size=15, bold=True, color=ACCENT)
    add_text(slide, "离焦参数 g = dz · rᵢ² / (2 · Ro · f)", top=1.9, size=13, color=GRAY)
    add_bullet(slide, [
        "g > 0 (焦点后): 形成焦散线 (caustic)，强度出现超越峰",
        "g < 0 (焦点前): 裙边宽度增加 Δa ≈ 3.88g − 3.53g²",
        "远离裙边区: 当 g ≪ 1 时影响可忽略",
        "裙边区内: 需要 β·g ≪ 1 才能保证影响可忽略",
    ], top=2.2, size=13)

    add_text(slide, "横向偏移 (d)：", top=4.1, size=15, bold=True, color=ACCENT)
    add_text(slide, "引入线性相位斜率，强度轮廓变为 |G(a)|² = 1 + 2·d·a", top=4.5, size=13, color=GRAY)
    add_bullet(slide, [
        "截断效应: 当入射光束被孔径截断时，裙边振荡幅度受限于 exp(-ξ₀²/2) · √(2β/π)",
        "此分析与本项目第二阶段的容差评估互补验证",
    ], top=4.9, size=13)

    # ====================================================================
    # Slide 6: RTAD — concept
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ② — RTAD 目标函数", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, "Chen et al., \"...high uniformity flat-top beams by RTAD,\" Optics & Laser Technology, 186, 112776 (2025)", top=0.85, size=10, color=GRAY)
    add_bullet(slide, [
        "问题: 传统平顶目标使用陡峭几何边界 → FFT 迭代产生频谱泄露 (Gibbs 现象) → 散斑噪声",
        "方案: 在平顶区与背景之间引入升余弦 (raised-cosine) 下降边缘",
    ], top=1.15, size=14)

    add_text(slide, "升余弦边缘函数 C(u)：", top=2.2, size=15, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"C(u)=1\ (u\leq u_0);\quad "
        r"C(u)=\frac{1}{2}\!\left[1+\cos\!\left(\pi\frac{u-u_0}{u_1-u_0}\right)\right]\ (u_0<u<u_1);\quad "
        r"C(u)=0\ (u\geq u_1)",
        left=0.5, top=2.6, fontsize=12, dpi=180)

    add_text(slide, "二维可分离目标强度:", top=3.7, size=14, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"I_{\mathrm{full}}(x,y) = C(|x|;\, a_0, a_1) \;\cdot\; C(|y|;\, b_0, b_1), \qquad "
        r"A_{\mathrm{full}} = \sqrt{I_{\mathrm{full}}}",
        left=0.6, top=4.05, fontsize=13, dpi=180)

    add_bullet(slide, [
        "a₀ = a₅₀ − Δx,  a₁ = a₅₀ + Δx    (默认 Δx=15μm, Δy=8μm)",
        "文献报告: 仿真均匀性 98.8%,  实验 91.0%,  显著优于传统几何目标法",
    ], top=4.8, size=13)

    # ====================================================================
    # Slide 7: RTAD — truncated constraint mode
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ② — 截断 RTAD 约束模式", top=0.3, size=28, bold=True, color=BLUE)

    add_text(slide, "传统 RTAD: 约束所有 I_full > 0 的像素 → 包含极低强度尾迹 → 过度约束", top=1.0, size=14, bold=True)
    add_text(slide, "本项目: 截断 RTAD (truncated RTAD) — 只约束强度足够高的区域", top=1.5, size=14, bold=True, color=ACCENT)

    add_bullet(slide, [
        f"release_level = exp(−2) ≈ 13.5%    (1/e² 强度水平)",
    ], top=2.0, size=14)

    add_text(slide, "焦平面区域划分:", top=2.5, size=15, bold=True, color=ACCENT)
    add_table(slide,
              ["区域", "定义", "投影方式"],
              [
                  ["mask_flat", "|x|≤a₀ ∩ |y|≤b₀   (纯平顶)", "WGS 权重反馈, 跟踪均匀性"],
                  ["mask_signal", "I_full ≥ 13.5%  OR  mask_flat", "强制目标振幅约束"],
                  ["mask_free", "保护窗口内非信号区", "MRAF 自由区 (mraf_factor·E)"],
                  ["mask_bg", "以上之外的全部", "远背景 (bg_factor·E, bg→1.0=不约束)"],
              ],
              col_widths=[1.8, 3.2, 3.8], font_size=11, top=3.0, row_height=0.42)

    add_bullet(slide, [
        "优点: 低强度尾迹进入自由区，不被强制约束 → 保留物理旁瓣/光晕，避免过度平滑",
        "mask_flat 是纯平顶区 (无边缘过渡)，RMS 均匀性仅在此区域计算",
    ], top=5.1, size=13)

    # ====================================================================
    # Slide 8: WGS — algorithm
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ③ — WGS 优化算法", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, "Alsaka et al., Applied Physics B, 128, 137 (2022)  |  算法参考: slmsuite v0.4.1", top=0.85, size=10, color=GRAY)

    add_text(slide, "迭代框架 (GS → MRAF → WGS):", top=1.1, size=15, bold=True, color=ACCENT)
    add_bullet(slide, [
        "①  前向 FFT: DOE 面复振幅 → 焦面复振幅",
        "②  输出面振幅投影 (MRAF): 信号区→目标振幅, 自由区→保留, 背景→弱衰减",
        "③  逆向 FFT: 焦面 → DOE 面, 提取新相位",
        "④  WGS 权重更新: 根据当前焦面振幅在 mask_flat 内自适应调整目标权重",
    ], top=1.5, size=13)

    add_text(slide, "MRAF 输出面投影:", top=3.2, size=14, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"E'_{\mathrm{signal}} = A_{\mathrm{target}} \cdot e^{i\cdot\arg(E)}, \quad "
        r"E'_{\mathrm{free}} = m_{\mathrm{mraf}} \cdot E, \quad "
        r"E'_{\mathrm{bg}} = b_{\mathrm{bg}} \cdot E",
        left=0.4, top=3.7, fontsize=14, dpi=180)

    add_text(slide, "WGS 权重反馈 (仅在 mask_flat 内更新):", top=4.3, size=14, bold=True, color=ACCENT)
    add_eq_image(slide,
        r"w \leftarrow w \cdot \left(\frac{\langle|E_{\mathrm{flat}}|\rangle}{|E|}\right)^{\gamma},\quad "
        r"w = \mathrm{clip}(w,\,w_{\min},\,w_{\max}),\quad "
        r"w\ /\!\!=\ \langle w \rangle",
        left=0.3, top=4.8, fontsize=14, dpi=180)

    add_bullet(slide, [
        "权重仅更新 mask_flat; 边缘、自由区、背景保持权重 1.0",
    ], top=5.6, size=12)

    # ====================================================================
    # Slide 9: WGS — key strategy + baseline
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "物理基础 ③ — 优化策略与基线结果", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, "最终优化流程: phase0 (Point2P) → direct flat-local WGS (不经 MRAF 预热)", top=0.85, size=14, bold=True, color=ACCENT)

    add_text(slide, "核心经验:", top=1.4, size=15, bold=True, color=BLUE)
    add_bullet(slide, [
        "★ 不约束远背景 (bg_factor→1.0) 是最大的单项改进 — 平台均匀性和剖面显著改善",
        "★ 直接 WGS 优于 MRAF-only (均匀性 80%→95%) 和 MRAF→WGS 混合",
        "★ release_level 在 5~13.5% 范围内影响很小; 反馈指数 0.5 以上饱和",
    ], top=1.8, size=13)

    add_text(slide, "最终工作基线参数:", top=3.2, size=15, bold=True, color=ACCENT)
    add_table(slide,
              ["参数", "值"],
              [
                  ["method", "wgs (不经过 MRAF)"],
                  ["strategy", "flat_local (2D 局部反馈)"],
                  ["iterations", "200"],
                  ["feedback_exponent", "0.8"],
                  ["weight range", "[0.5, 2.0]"],
                  ["mraf_factor (自由区)", "0.8"],
                  ["bg_factor (远背景)", "1.0 (实质不约束)"],
                  ["release_level", "0.135 (exp(-2))"],
              ],
              col_widths=[4.4, 4.4], font_size=12, top=3.8, row_height=0.32)

    add_text(slide, "RMS 非均匀性 = 1.86%   |   e⁻² 衍射效率 = 92.5%   |   size50 = 330.2 × 123.6 μm",
             top=6.7, size=13, bold=True, color=ACCENT)

    # ====================================================================
    # Slide 10: Tolerance assessment
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "真实光路容差评估 — 固定 DOE 相位，模拟误差感度", top=0.3, size=26, bold=True, color=BLUE)
    add_text(slide, "目标: 评估已设计好的 DOE 在真实光路条件下的鲁棒性，指导实验反向纠错", top=0.85, size=14, bold=True, color=ACCENT)

    add_table(slide,
              ["误差类型", "参数", "物理含义"],
              [
                  ["离焦 (defocus)", "defocus_mm", "观察面偏离焦平面 ±2mm"],
                  ["光束偏移 (offset)", "offset_x/y_mm", "光斑偏离 DOE 中心 ±0.5mm"],
                  ["光束尺寸 (beam size)", "diameter_1e2_mm", "入射 1/e² 直径偏差 3.5~6.5mm"],
                  ["发散/会聚 (divergence)", "divergence_edge_mrad", "DOE 面波前曲率 ±0.5mrad"],
                  ["倾角 (pointing)", "pointing_shift_x/y_um", "整束倾斜 → 焦斑平移 ±100μm"],
                  ["通光孔径 (aperture)", "clear_aperture_mm", "有效通光区域 10~15mm"],
                  ["椭圆度 (ellipticity)", "Dx / Dy", "x/y 光斑尺寸不一致"],
              ],
              col_widths=[2.5, 2.5, 3.8], font_size=13, top=1.7, row_height=0.60)

    add_text(slide, "每种误差独立扫描 (mild + stress)，双参数组合扫描用于定位耦合效应", top=6.5, size=13, color=GRAY)

    # ====================================================================
    # Slide 11: Error fitting results
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "5 种典型光斑误差 — 根因拟合结果", top=0.3, size=26, bold=True, color=BLUE)
    add_table(slide,
              ["#", "问题现象", "根因", "物理机制"],
              [
                  ["①", "长边中间内凹", "椭圆度 Dx > Dy", "x 方向过度照明 → 抢走长边中心处 y 向能量"],
                  ["②", "长边能量缺失 + 短边倾斜", "beam_offset (x+y 复合偏移)", "光斑偏离中心造成不对称照明"],
                  ["③", "四条边都向内凹", "负离焦 + 光斑偏大 (≥5.5mm)", "观察面在焦点前 + 扩束过大 → 四边对称内缩"],
                  ["④", "能量上下聚集", "beam_size (光束尺寸偏差)", "入射光斑尺寸改变有效照明轮廓"],
                  ["⑤", "长宽比不对", "待进一步分析", "可能与椭圆度/发散角组合相关"],
              ],
              col_widths=[0.4, 2.5, 2.5, 3.4], font_size=13, top=1.3, row_height=0.90)

    add_text(slide, "容差量级参考 (RMS<5%) — 发散角 ~±0.005mrad | 离焦 ~±0.75mm | 偏移 ~±0.2mm | 椭圆度 <0.2mm",
             top=6.5, size=12, color=GRAY)

    # ====================================================================
    # Slide 12: Project structure
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "项目架构", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, """Point2P/
├── initial_phase_generation/    [MATLAB]  Romero-Dickey 初始相位生成 (已冻结)
│   └── generate_initial_phase.m    核心: 稳相法解析相位, 输出 phase0.mat
│
├── rtad_mraf_gs_python/         [Python]  MRAF/WGS 迭代优化
│   ├── src/mraf_gs.py              核心: GS/MRAF/WGS 精化循环
│   ├── src/rtad_target.py          RTAD 升余弦矩形目标构建
│   ├── src/propagation.py          FFT 传播 + 高斯入射场生成
│   └── run_rtad_mraf_gs_case.py    CLI 入口
│
├── real_world_simulation/       [Python]  真实光路容差评估
│   ├── src/field_models.py         真实入射场建模 (偏移/发散/倾角等)
│   ├── src/propagation.py          传播 + 角谱法离焦
│   ├── src/metrics.py              诊断指标 (size50/13.5/90, RMS, e²)
│   └── run_real_world_sweep.py     CLI 扫描入口
│
├── result_diagnostics/          [MATLAB]  焦平面诊断 (备用)
├── target/                      [MATLAB]  RTAD 目标生成 (备用)
└── text/                        参考文献 (×4)""", size=11, top=1.0, left=0.4, width=9.2, height=6.2)

    # ====================================================================
    # Slide 13: Summary
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "总结", top=0.3, size=28, bold=True, color=BLUE)
    add_bullet(slide, [
        "Point2P 初始相位: Romero-Dickey 解析解提供接近目标能量分布的物理初值 (βₓ=9.06, βⵧ=3.30)",
        "RTAD 目标函数: 升余弦下降边缘有效抑制 FFT 频谱泄露; 截断约束保留物理旁瓣结构",
        "WGS 优化: 直接 flat-local WGS, 不约束远背景 (bg_factor→1.0) 是最大单项改进",
        "最终结果: RMS 非均匀性 = 1.86%,  e⁻² 衍射效率 = 92.5%",
        "容差分析: 揭示 5 类光斑误差的根因，可用于实验反向纠错 (光斑形态 → 误差参数)",
    ], top=1.3, size=16)

    add_text(slide, "关键容差 — 发散角 ~±0.005 mrad (极敏感)  |  离焦 ~±0.75 mm  |  偏移 ~±0.2 mm", top=6.5, size=13, color=ACCENT)

    # ====================================================================
    # Slide 14: References
    # ====================================================================
    slide = prs.slides.add_slide(slide_layout(prs, 5))
    add_text(slide, "参考文献", top=0.3, size=28, bold=True, color=BLUE)
    add_text(slide, """[1]  Romero, L.A. & Dickey, F.M.
     "Lossless laser beam shaping."  J. Opt. Soc. Am. A, 13(4), 751–760 (1996).
     → 初始相位理论基础: 稳相法解析求解高斯→平顶 DOE 相位

[2]  Chen, W. et al.
     "Generation of high uniformity flat-top beams by reconstructing
      the amplitude distribution at descending edges."
     Optics & Laser Technology, 186, 112776 (2025).
     → RTAD 目标函数: 升余弦下降边缘替代陡峭几何边界

[3]  Alsaka, D.Y. et al.
     "Dynamic flat-topped laser beam shaping method using mixed region
      amplitude freedom algorithm."  Applied Physics B, 128, 137 (2022).
     → MRAF/WGS 算法: 混合区域振幅自由 + 自适应权重反馈

[4]  Zhang, C. et al.
     "Optimized holographic femtosecond laser patterning method towards
      rapid integration of high-quality functional devices in microchannels."
     Scientific Reports, 6, 33281 (2016).
     → MRAF 在飞秒激光微加工中的应用

[5]  slmsuite v0.4.1  —  github.com/slmsuite/slmsuite
     → 算法实现参考: FFT 传播、振幅归一化、MRAF/WGS 语义""", size=12, top=1.2, left=0.5, width=9.0, height=6.0)

    # ====================================================================
    prs.save(OUTPUT)
    print(f"Done: {OUTPUT}")


if __name__ == "__main__":
    make_ppt()
