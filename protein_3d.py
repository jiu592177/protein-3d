import streamlit as st
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="3D 蛋白质合成", layout="wide")
st.title("🧬 3D 蛋白质合成演示")

# ================= 数据 =================
dna = "TACGGCTTAACGCAT"
comp = {"A":"U", "T":"A", "C":"G", "G":"C"}
mrna_full = "".join(comp[b] for b in dna)
n = len(dna)

dna_colors = {"A":"#e74c3c", "T":"#3498db", "C":"#27ae60", "G":"#f39c12"}
mrna_colors = {"A":"#e74c3c", "U":"#9b59b6", "C":"#27ae60", "G":"#f39c12"}

codon_table = {
    "AUG":"Met", "CCG":"Pro", "AAU":"Asn", "UGC":"Cys", "GUA":"Val",
    "UAA":"Stop", "UAG":"Stop", "UGA":"Stop",
}
codons = [mrna_full[i:i+3] for i in range(0, n-2, 3)]
aas = [codon_table.get(c, "?") for c in codons]

# ================= 公共函数 =================
def helix_coords(n_bases, radius=1.0, rise=0.75):
    turns = 0.35
    theta = np.linspace(0, 2*np.pi*turns*n_bases, n_bases*10)
    z = np.linspace(0, rise*n_bases, n_bases*10)
    x1, y1 = radius*np.cos(theta), radius*np.sin(theta)
    x2, y2 = radius*np.cos(theta+np.pi), radius*np.sin(theta+np.pi)
    return x1, y1, z, x2, y2

def base_layout(fig, title=""):
    fig.update_layout(
        scene=dict(
            xaxis=dict(showgrid=False, showticklabels=False, backgroundcolor="#f8f9fa"),
            yaxis=dict(showgrid=False, showticklabels=False, backgroundcolor="#f8f9fa"),
            zaxis=dict(title="高度", showgrid=True, gridcolor="#ecf0f1",
                       backgroundcolor="#f8f9fa"),
            camera=dict(eye=dict(x=1.8, y=-1.8, z=1.0)),
            aspectmode="manual",
            aspectratio=dict(x=1.2, y=1, z=1.5),
        ),
        legend=dict(x=0, y=1, bgcolor="rgba(255,255,255,0.8)"),
        height=680,
        margin=dict(l=0, r=0, t=30, b=0),
        paper_bgcolor="white",
        title=dict(text=title, x=0.5, font=dict(size=16)),
    )
    return fig

# ================= 建 Tab 分页 =================
tab1, tab2, tab3 = st.tabs(["① 转录过程", "② 翻译过程", "③ 蛋白质折叠"])

# =====================================================
#                  Tab 1：转录
# =====================================================
with tab1:
    col_ctrl, col_plot = st.columns([1, 4])

    with col_ctrl:
        st.subheader("控制")
        transcribe_pos = st.slider("转录进度", 0, n, 0, 1, key="t_slider")
        st.markdown("---")
        st.caption("RNA 聚合酶沿 DNA 模板链移动，\nmRNA 逐个碱基延伸。")
        st.metric("已转录", f"{transcribe_pos} / {n}")

    with col_plot:
        fig = go.Figure()
        x1, y1, z, x2, y2 = helix_coords(n)

        # DNA
        fig.add_trace(go.Scatter3d(x=x1, y=y1, z=z, mode="lines",
            line=dict(color="#3498db", width=8), name="DNA 模板链", opacity=0.85))
        fig.add_trace(go.Scatter3d(x=x2, y=y2, z=z, mode="lines",
            line=dict(color="#2980b9", width=8), name="DNA 编码链", opacity=0.85))

        for i in range(n):
            idx = i * 10
            if idx < len(x1):
                fig.add_trace(go.Scatter3d(
                    x=[x1[idx], x2[idx]], y=[y1[idx], y2[idx]], z=[z[idx], z[idx]],
                    mode="lines", line=dict(color="#95a5a6", width=3),
                    showlegend=False, opacity=0.6))
                fig.add_trace(go.Scatter3d(
                    x=[x1[idx]], y=[y1[idx]], z=[z[idx]],
                    mode="markers+text",
                    marker=dict(size=7, color=dna_colors[dna[i]]),
                    text=[dna[i]], textposition="middle right",
                    textfont=dict(size=9, color=dna_colors[dna[i]]),
                    showlegend=False, hoverinfo="skip"))

        # mRNA
        if transcribe_pos > 0:
            mx, my, mz = [], [], []
            for i in range(transcribe_pos):
                idx = i * 10
                if idx < len(x1):
                    mx.append(x1[idx] + 0.6)
                    my.append(y1[idx] + 0.6)
                    mz.append(z[idx] - 0.5)
            if mx:
                fig.add_trace(go.Scatter3d(
                    x=mx, y=my, z=mz,
                    mode="lines+markers+text",
                    line=dict(color="#e74c3c", width=7),
                    marker=dict(size=9,
                        color=[mrna_colors[mrna_full[i]] for i in range(len(mx))],
                        line=dict(color="white", width=1)),
                    text=[mrna_full[i] for i in range(len(mx))],
                    textposition="bottom center",
                    textfont=dict(size=11, color="#c0392b"),
                    name="mRNA"))

                # 高亮最新碱基
                if transcribe_pos < n:
                    fig.add_trace(go.Scatter3d(
                        x=[mx[-1]], y=[my[-1]], z=[mz[-1]],
                        mode="markers",
                        marker=dict(size=18, color="#e74c3c", opacity=0.4),
                        showlegend=False, hoverinfo="skip"))

        # RNA 聚合酶
        if transcribe_pos < n:
            idx = transcribe_pos * 10
            if idx < len(x1):
                px, py, pz = x1[idx], y1[idx], z[idx]
                fig.add_trace(go.Mesh3d(
                    x=px + 0.5*np.cos(np.linspace(0, 2*np.pi, 25)),
                    y=py + 0.5*np.sin(np.linspace(0, 2*np.pi, 25)),
                    z=pz + 0.5*np.sin(np.linspace(0, 2*np.pi, 25)),
                    color="#8e44ad", opacity=0.65, name="RNA 聚合酶"))

        base_layout(fig, "转录：DNA → mRNA")
        st.plotly_chart(fig, use_container_width=True)

    st.markdown("### 🧬 序列对照")
    st.code(
        f"DNA 模板链 : 5'-{dna}-3'\n"
        f"mRNA       : 5'-{mrna_full[:transcribe_pos]}{'·'*(n-transcribe_pos)}-3'"
    )

# =====================================================
#                  Tab 2：翻译
# =====================================================
with tab2:
    col_ctrl, col_plot = st.columns([1, 4])

    with col_ctrl:
        st.subheader("控制")
        translate_pos = st.slider("翻译进度（氨基酸数）",
                                  0, len(aas), 0, 1, key="tr_slider")
        st.markdown("---")
        st.caption("核糖体读取密码子，\ntRNA 运送氨基酸，形成肽键。")
        st.metric("已翻译", f"{translate_pos} / {len(aas)}")

    with col_plot:
        fig = go.Figure()
        x1, y1, z, x2, y2 = helix_coords(n)

        # DNA（暗一点当背景）
        fig.add_trace(go.Scatter3d(x=x1, y=y1, z=z, mode="lines",
            line=dict(color="#3498db", width=5), name="DNA", opacity=0.3))
        fig.add_trace(go.Scatter3d(x=x2, y=y2, z=z, mode="lines",
            line=dict(color="#2980b9", width=5), name="DNA", opacity=0.3,
            showlegend=False))

        # mRNA 完整
        mx, my, mz = [], [], []
        for i in range(n):
            idx = i * 10
            if idx < len(x1):
                mx.append(x1[idx] + 0.6)
                my.append(y1[idx] + 0.6)
                mz.append(z[idx] - 0.5)

        fig.add_trace(go.Scatter3d(
            x=mx, y=my, z=mz,
            mode="lines+markers+text",
            line=dict(color="#e74c3c", width=7),
            marker=dict(size=9,
                color=[mrna_colors[mrna_full[i]] for i in range(len(mx))],
                line=dict(color="white", width=1)),
            text=list(mrna_full),
            textposition="bottom center",
            textfont=dict(size=11, color="#c0392b"),
            name="mRNA"))

        # 氨基酸链
        n_aa = int(translate_pos)
        if n_aa > 0:
            aa_x = mx[:n_aa*3:3][:n_aa]
            aa_y = np.full(n_aa, -1.5)
            aa_z = np.full(n_aa, 1.5)

            if n_aa > 1:
                fig.add_trace(go.Scatter3d(
                    x=aa_x, y=aa_y, z=aa_z,
                    mode="lines", line=dict(color="#7f8c8d", width=4),
                    name="肽键", showlegend=False))

            colors = ["#f39c12","#e67e22","#f1c40f","#d35400",
                      "#e74c3c","#9b59b6","#3498db","#1abc9c"]
            fig.add_trace(go.Scatter3d(
                x=aa_x, y=aa_y, z=aa_z,
                mode="markers+text",
                marker=dict(size=12,
                    color=[colors[i%len(colors)] for i in range(n_aa)]),
                text=aas[:n_aa], textposition="bottom center",
                textfont=dict(size=10, color="#7d3c98"),
                customdata=codons[:n_aa],
                hovertemplate="密码子: %{customdata}<br>氨基酸: %{text}<extra></extra>",
                name="氨基酸链"))

            # 核糖体
            li = min(n_aa*3, len(mx)-1)
            fig.add_trace(go.Mesh3d(
                x=mx[li] + 0.45*np.cos(np.linspace(0, 2*np.pi, 25)),
                y=0.45*np.sin(np.linspace(0, 2*np.pi, 25)),
                z=mz[li] + 0.8 + 0.45*np.sin(np.linspace(0, 2*np.pi, 25)),
                color="#27ae60", opacity=0.6, name="核糖体"))

        base_layout(fig, "翻译：mRNA → 多肽链")
        st.plotly_chart(fig, use_container_width=True)

    st.markdown("### 🧬 序列对照")
    st.code(
        f"mRNA   : 5'-{mrna_full}-3'\n"
        f"多肽链 : {'-'.join(aas[:int(translate_pos)])}"
    )

# =====================================================
#                  Tab 3：折叠（独立页面）
# =====================================================
with tab3:
    col_ctrl, col_plot = st.columns([1, 4])

    with col_ctrl:
        st.subheader("控制")
        fold_pos = st.slider("折叠进度", 0, 100, 0, 1, key="f_slider")
        st.markdown("---")
        st.caption("多肽链从伸展状态\n折叠成球状三维结构。")
        st.metric("折叠进度", f"{fold_pos}%")

        st.markdown("**结构层次**")
        st.markdown("- 一级：氨基酸排列顺序\n"
                    "- 二级：α螺旋 / β折叠\n"
                    "- 三级：整条链的空间构象\n"
                    "- 四级：多条链组装")

    with col_plot:
        fig = go.Figure()
        n_aa = len(aas)
        t = fold_pos / 100.0

        # 生成折叠后的球面点
        golden = np.pi * (3 - np.sqrt(5))
        sphere_pts = []
        for i in range(n_aa):
            y_s = 1 - (i / max(n_aa - 1, 1)) * 2
            r_s = np.sqrt(max(1 - y_s*y_s, 0))
            theta_s = golden * i
            sphere_pts.append([np.cos(theta_s)*r_s, y_s, np.sin(theta_s)*r_s])
        sphere_pts = np.array(sphere_pts)

        # 线性伸展位置（横排）
        linear_pos = np.column_stack([
            np.arange(n_aa) * 1.2 - (n_aa-1)*0.6,
            np.zeros(n_aa),
            np.zeros(n_aa)
        ])

        # 折叠后位置（球形）
        center = np.array([0, 0, 0])
        folded_pos = center + sphere_pts * 2.0

        # 插值
        cur_pos = (1 - t) * linear_pos + t * folded_pos

        colors = ["#f39c12","#e67e22","#f1c40f","#d35400",
                  "#e74c3c","#9b59b6","#3498db","#1abc9c"]

        # 肽链骨架
        fig.add_trace(go.Scatter3d(
            x=cur_pos[:,0], y=cur_pos[:,1], z=cur_pos[:,2],
            mode="lines",
            line=dict(color="#95a5a6", width=4),
            name="肽链骨架", showlegend=False))

        # 氨基酸
        fig.add_trace(go.Scatter3d(
            x=cur_pos[:,0], y=cur_pos[:,1], z=cur_pos[:,2],
            mode="markers+text",
            marker=dict(size=14,
                color=[colors[i%len(colors)] for i in range(n_aa)],
                line=dict(color="white", width=1)),
            text=aas, textposition="bottom center",
            textfont=dict(size=10, color="#7d3c98"),
            name="氨基酸",
            hovertemplate="氨基酸: %{text}<extra></extra>"))

        # 折叠完成后的“蛋白质外壳”
        if t > 0.85:
            u = np.linspace(0, 2*np.pi, 40)
            v = np.linspace(0, np.pi, 20)
            U, V = np.meshgrid(u, v)
            Xs = 2.6 * np.cos(U) * np.sin(V)
            Ys = 2.6 * np.sin(U) * np.sin(V)
            Zs = 2.6 * np.cos(V)
            fig.add_trace(go.Surface(
                x=Xs, y=Ys, z=Zs,
                colorscale=[[0, "rgba(39,174,96,0.0)"],
                            [1, "rgba(39,174,96,0.25)"]],
                showscale=False, name="蛋白质外壳",
                hoverinfo="skip"))

        # 状态文字
        label = "✅ 折叠完成，具有生物活性的蛋白质" if t > 0.85 \
                else ("伸展的多肽链" if t < 0.15 else f"折叠中 {fold_pos}%")
        fig.add_trace(go.Scatter3d(
            x=[0], y=[0], z=[3.8],
            mode="text", text=[label],
            textfont=dict(size=15, color="#27ae60"),
            showlegend=False))

        fig.update_layout(
            scene=dict(
                xaxis=dict(showgrid=False, showticklabels=False, backgroundcolor="#f8f9fa"),
                yaxis=dict(showgrid=False, showticklabels=False, backgroundcolor="#f8f9fa"),
                zaxis=dict(showgrid=False, showticklabels=False, backgroundcolor="#f8f9fa"),
                camera=dict(eye=dict(x=1.5, y=-1.5, z=1.2)),
                aspectmode="cube",
            ),
            height=680,
            margin=dict(l=0, r=0, t=30, b=0),
            paper_bgcolor="white",
            title=dict(text="蛋白质折叠：多肽链 → 三维结构",
                       x=0.5, font=dict(size=16)),
        )
        st.plotly_chart(fig, use_container_width=True)

    st.markdown("### 🧬 折叠说明")
    st.code(
        f"多肽链 : {'-'.join(aas)}\n"
        f"长度   : {n_aa} 个氨基酸"
    )
    st.info("拖动左侧滑块，观察多肽链从线性伸展逐渐蜷曲成球状蛋白质。")