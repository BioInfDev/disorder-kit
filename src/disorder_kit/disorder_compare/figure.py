import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

def bar_disorder_distribution(bins, bin_colors, method):
    counts = bins.value_counts().sort_index()
    labels = counts.index.astype(str).tolist()
    values = counts.values.tolist()

    fig = go.Figure(go.Bar(
            x = labels,
            y = values,
            text = values,
            textposition = 'outside',
            textfont = dict(size = 14, color = '#333', weight = 'bold'),
            marker = dict(color = bin_colors),
            width = 0.7,
        )
    )
    fig.update_traces(marker_cornerradius=25)
    fig.update_layout(
        title=dict(text = "<b>Distribution of proteins by PPIDR (Percentage of Protein Intrinsic Disorder Residues)</b>"
                        f"<br><span style='font-size:12px;color:gray'>Disorder was predicted using {method}</span>",
                    x = 0.1),
        yaxis = dict(visible = False),
        xaxis = dict(showgrid = False, zeroline = False, tickfont = dict(size=12, color='#666')),
        plot_bgcolor = 'white',
        showlegend = False,
        width = 1000
    )
    return fig

def scatter_disorder_distribution(df, x_col, y_col, bins, bin_colors, bins_label, targets, hover_cols):
    fig = px.scatter(
        df, 
        x = x_col, 
        y = y_col,
        color = bins,
        color_discrete_sequence = bin_colors,
        labels = {'color': bins_label},
        category_orders={'color': bins.cat.categories},
        opacity = 0.8,
        hover_data = hover_cols,
        template = 'plotly_white',
        title = "<b>Distribution of proteins by intrinsic disorder content</b><br>"
            "<span style='font-size:12px;color:gray'>"
            "Disorder was predicted using PPIDR-AIUPred and PPIDR-pLDDT</span>"
    )

    fig.update_traces(
        marker=dict(size=8, line=dict(width=0)),
        hovertemplate=f"<b>{hover_cols[0]}</b><br>"
                    f"{hover_cols[1]}: %{{x:.1f}}%<br>"
                    f"{hover_cols[2]}: %{{y:.1f}}%<br>"
                    f"{hover_cols[3]}: %{{customdata[1]:.1f}}"
                    f"{hover_cols[4]}: %{{customdata[2]:.1f}}"
                    "<extra></extra>"
    )

    color_map = {tr.name: tr.marker.color for tr in fig.data if tr.name}
    target_bins = bins.loc[targets.index]
    marker_colors = target_bins.map(color_map).fillna('#888888').tolist()

    fig.add_trace(go.Scatter(
        x = targets[x_col],
        y = targets[y_col],
        mode = 'markers+text',
        marker = dict(
            size = 16,
            color = marker_colors,
            line = dict(width = 2, color = 'white'),
            symbol = 'circle',
        ),
        text=targets[df.columns[0]],
        textposition = 'top center',
        textfont = dict(size = 15, color = '#333'),
        showlegend = False,                     
        customdata = targets[[df.columns[0], df.columns[3]]].values,
        hovertemplate=f"<b>%{hover_cols[0]}</b><br>"
                    f"{hover_cols[1]}: %{{x:.1f}}%<br>"
                    f"{hover_cols[2]}: %{{y:.1f}}%<br>"
                    f"{hover_cols[3]}: %{{customdata[3]:.1f}}"
                    f"{hover_cols[4]}: %{{customdata[4]:.1f}}"
                    "<extra></extra>"
    ))
    return fig

def sankey_distribution_compare(compare_bins, titles):   
    bin_col1 = compare_bins[0]
    bin_col2 = compare_bins[1]
    
    labels = bin_col1.cat.categories
    n = len(labels)
    
    if titles is None:
        titles = ["Input 1", "Input 2"]
        
    # Маппинги для двух колонок (смещаем индексы правой стороны на n)
    map1 = {val: i for i, val in enumerate(labels)}
    map2 = {val: i + n for i, val in enumerate(labels)}
    
    all_labels = (
        [f"{titles[0]}: {l}" for l in labels] + 
        [f"{titles[1]}: {l}" for l in labels]
    )
    
    color_left = "#A8BEDC"   # Левая сторона
    color_right = "#409151"  # Правая сторона
    
    node_colors = [color_left] * n + [color_right] * n
    
    df_12 = pd.DataFrame({'c1': bin_col1, 'c2': bin_col2})
    counts_12 = df_12.groupby(['c1', 'c2'], observed=False).size().reset_index(name='count')
    counts_12 = counts_12[counts_12['count'] > 0]
    
    sources = counts_12['c1'].map(map1).tolist()
    targets = counts_12['c2'].map(map2).tolist()
    values = counts_12['count'].tolist()
    
    link_color = "rgba(168, 190, 220, 0.3)"
    link_colors = [link_color] * len(sources)
    
    # Отрисовка графика
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="rgba(0,0,0,0.3)", width=0.5),
            label=all_labels,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors
        )
    )])

    fig.update_layout(
        title = dict(
            text = f"Overlap beetween '{titles[0]}' and '{titles[1]}'",
            x = 0.1),
        font_size=13,
        width=1500,
        height=700
    )
    return fig