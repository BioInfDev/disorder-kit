import plotly.graph_objects as go
import plotly.express as px


def bar_disorder_distribution(bins, title):
    counts = bins.value_counts().sort_index()
    labels = counts.index.astype(str).tolist()
    values = counts.values.tolist()

    fig = go.Figure(go.Bar(
        x = labels,
        y = values,
        text = values,
        textposition = 'outside',
        textfont=dict(size = 14, color = '#333', weight = 'bold'),
        marker = dict(color=['#5A5DF6', '#8E93F9', '#B0B0B0', '#C9C9C9', '#E0E0E0']),
        width = 0.5
    ))

    fig.update_traces(marker_cornerradius=15)
    fig.update_layout(
        title=dict(
            text = "<b>Distribution of proteins by intrinsic disorder content</b>"
                f"<br><span style='font-size:12px;color:gray'>Disorder was predicted using {title}</span>",
            x = 0
        ),
        yaxis = dict(showgrid=False, rangemode='tozero'),
        xaxis = dict(showgrid=False, zeroline=False,
                tickfont=dict(size=12, color='#666')),
        plot_bgcolor = 'white',
        showlegend = False,
        margin = dict(t=80, b=40, l=40, r=20)
    )
    fig.show()

def scatter_disorder_distribution(df, df_targets, chm_bins, legend_label: str = None):
    fig = px.scatter(
        df, 
        x = df.columns[1], 
        y = df.columns[2],
        color = chm_bins,
        category_orders = {'color': ['0–20', '20–40', '40–60', '60–80', '80–100']},
        labels = {'color': df.columns[3]},
        color_continuous_scale = 'Viridis',
        opacity = 0.7,
        hover_data = [df.columns[0], df.columns[3]],
        template = 'plotly_white',
        title = "<b>Distribution of proteins by intrinsic disorder content</b><br>"
            "<span style='font-size:12px;color:gray'>"
            "Disorder was predicted using AIUPred and pLDDT</span>"
    )

    fig.update_traces(
        marker=dict(size=8, line=dict(width=0)),
        hovertemplate="<b>%{customdata[0]}</b><br>"
                    f"{df.columns[1]}: %{{x:.1f}}%<br>"
                    f"{df.columns[2]}: %{{y:.1f}}%<br>"
                    f"{df.columns[3]}: %{{customdata[1]:.1f}}"
                    "<extra></extra>",
    )

    color_map = {tr.name: tr.marker.color for tr in fig.data if tr.name}
    target_bins = chm_bins.loc[df_targets.index]
    marker_colors = target_bins.map(color_map).fillna('#888888').tolist()

    fig.add_trace(go.Scatter(
        x = df_targets[df.columns[1]],
        y = df_targets[df.columns[2]],
        mode = 'markers+text',
        marker = dict(
            size = 16,
            color = marker_colors,
            line = dict(width = 2, color = 'white'),
            symbol = 'circle',
        ),
        text=df_targets[df.columns[0]],
        textposition = 'top center',
        textfont = dict(size = 15, color = '#333'),
        showlegend = False,                     
        customdata = df_targets[[df.columns[0], df.columns[3]]].values,
        hovertemplate="<b>%{customdata[0]}</b><br>"
                            f"{df.columns[1]}: %{{x:.1f}}%<br>"
                            f"{df.columns[2]}: %{{y:.1f}}%<br>"
                            f"{df.columns[3]}: %{{customdata[1]:.1f}}"
                            "<extra></extra>",
    ))
    fig.show()