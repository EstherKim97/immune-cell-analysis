"""Streamlit dashboard for the Loblaw Bio assessment."""

import math
import sqlite3
from contextlib import closing

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

from frequency_table import DB_PATH, get_frequency_table
from response_analysis import run_analysis
from subset_analysis import (get_baseline_cohort,
                             get_baseline_male_responder_b_cell_average,
                             get_summary_tables)


POPULATION_LABELS = {
    "b_cell": "B cells",
    "cd8_t_cell": "CD8+ T cells",
    "cd4_t_cell": "CD4+ T cells",
    "nk_cell": "NK cells",
    "monocyte": "Monocytes",
}
POPULATION_ORDER = tuple(POPULATION_LABELS)
CRIMSON = "#8F263D"
DARK_CRIMSON = "#65192A"
PALE_CRIMSON = "#F7ECEF"
CHARCOAL = "#202124"
MUTED = "#6E6A69"
BORDER = "#DED9D6"
GRID = "#EAE6E3"
GROUP_COLORS = {"no": "#92908F", "yes": CRIMSON}
PAGES = ("01  Study Overview", "02  Immune Composition",
         "03  Treatment Response", "04  Cohort Explorer")


def connect():
    return sqlite3.connect(DB_PATH.as_uri() + "?mode=ro", uri=True)


@st.cache_data(show_spinner=False)
def load_overview(database_revision: tuple[int, int]) -> dict:
    with closing(connect()) as connection:
        counts = connection.execute("""
            SELECT (SELECT COUNT(*) FROM projects),
                   (SELECT COUNT(*) FROM subjects),
                   (SELECT COUNT(*) FROM samples),
                   (SELECT COUNT(*) FROM cell_counts),
                   (SELECT COUNT(DISTINCT cell_type) FROM cell_counts)
        """).fetchone()
        sample_counts = connection.execute("""
            SELECT sample_count, COUNT(*) FROM (
                SELECT COUNT(*) AS sample_count FROM samples GROUP BY subject_id
            ) GROUP BY sample_count ORDER BY sample_count
        """).fetchall()
    return dict(zip(("projects", "subjects", "samples", "measurements", "populations"),
                    counts)) | {"samples_per_subject": sample_counts}


@st.cache_data(show_spinner="Loading sample frequencies…")
def load_frequencies(database_revision: tuple[int, int]) -> pd.DataFrame:
    with closing(connect()) as connection:
        return pd.DataFrame.from_records(
            get_frequency_table(connection),
            columns=("sample", "total_count", "population", "count", "percentage"),
        )


@st.cache_data(show_spinner="Loading response analysis…")
def load_response_analysis(database_revision: tuple[int, int]) -> dict:
    with closing(connect()) as connection:
        return run_analysis(connection)


@st.cache_data(show_spinner=False)
def load_baseline_subset(database_revision: tuple[int, int]) -> dict:
    with closing(connect()) as connection:
        cohort = get_baseline_cohort(connection)
        return {
            "cohort": cohort,
            "tables": get_summary_tables(cohort),
            "b_cell_query": get_baseline_male_responder_b_cell_average(connection),
        }


def page_header(number: str, title: str, question: str) -> None:
    st.markdown(f'<div class="eyebrow">ANALYSIS {number}</div>', unsafe_allow_html=True)
    st.title(title)
    st.markdown(f'<p class="page-question">{question}</p><div class="page-rule"></div>',
                unsafe_allow_html=True)


def section_title(title: str, detail: str = "") -> None:
    st.markdown(f'<div class="section-heading"><h2>{title}</h2><span>{detail}</span></div>',
                unsafe_allow_html=True)


def metrics(items: list[tuple[str, str]]) -> None:
    columns = st.columns(len(items), gap="medium")
    for column, (label, value) in zip(columns, items):
        column.markdown(
            f'<div class="report-metric"><div class="metric-value">{value}</div>'
            f'<div class="metric-label">{label}</div></div>', unsafe_allow_html=True)


def chart_layout(figure: go.Figure, height: int, margin: dict | None = None) -> None:
    figure.update_layout(
        height=height, margin=margin or dict(l=40, r=18, t=22, b=30),
        paper_bgcolor="white", plot_bgcolor="white", font=dict(color=CHARCOAL, size=12),
        hoverlabel=dict(bgcolor="white", font_color=CHARCOAL),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    figure.update_yaxes(showgrid=True, gridcolor=GRID, zeroline=False,
                        linecolor=BORDER, ticks="outside", tickcolor=BORDER)
    figure.update_xaxes(showgrid=False, zeroline=False, linecolor=BORDER,
                        ticks="outside", tickcolor=BORDER)


def evidence_text(row) -> str:
    if row.fdr_q_value < .05:
        return "The difference remained statistically significant after FDR correction."
    if row.raw_p_value < .05:
        return "A nominal difference was observed, but it did not meet the FDR-adjusted significance threshold."
    return "No nominal evidence of a difference was observed."


def strongest(table: pd.DataFrame):
    return table.loc[table["raw_p_value"].idxmin()]


def overview_page(database_revision: tuple[int, int]) -> None:
    page_header("01", "Study Overview", "Clinical trial immune profiling and treatment-response analysis")
    st.markdown("Evaluate immune-cell composition across clinical trial samples and assess whether "
                "measured populations differ between melanoma patients who respond to miraclib "
                "and those who do not.")
    counts = load_overview(database_revision)
    section_title("Dataset at a glance")
    metrics([("Subjects", f"{counts['subjects']:,}"), ("Samples", f"{counts['samples']:,}"),
             ("Projects", f"{counts['projects']:,}"),
             ("Immune populations", f"{counts['populations']:,}")])
    section_title("Study structure", "Relational hierarchy")
    structure = st.columns(4, gap="small")
    for column, label, value in zip(structure,
                                    ("Projects", "Subjects", "Samples", "Cell measurements"),
                                    (counts["projects"], counts["subjects"], counts["samples"],
                                     counts["measurements"])):
        column.markdown(f'<div class="structure-step"><strong>{value:,}</strong>'
                        f'<span>{label}</span></div>', unsafe_allow_html=True)
    distribution = counts["samples_per_subject"]
    sample_note = (f"Every subject contributes {distribution[0][0]} longitudinal samples. "
                   if len(distribution) == 1 else
                   "Subjects contribute varying numbers of samples. ")
    st.caption(sample_note + "Sample-level measurements and subject-level counts are treated separately.")
    section_title("Analysis questions")
    for number, question in (
        ("01", "What fraction of each sample is represented by each immune-cell population?"),
        ("02", "Does immune-cell composition differ between miraclib responders and nonresponders?"),
        ("03", "What characterizes the baseline melanoma/miraclib/PBMC cohort?"),
    ):
        st.markdown(f'<div class="question-row"><span>{number}</span><p>{question}</p></div>',
                    unsafe_allow_html=True)
    section_title("Key result", "Pooled melanoma · miraclib · PBMC comparison")
    lead = strongest(load_response_analysis(database_revision)["pooled"])
    st.markdown(f'**Strongest observed difference: {POPULATION_LABELS[lead.population]}**')
    st.markdown(evidence_text(lead))
    st.caption(f"Responder median {lead.responder_median:.2f}% · nonresponder median "
               f"{lead.nonresponder_median:.2f}% · raw p = {lead.raw_p_value:.4g} · "
               f"FDR q = {lead.fdr_q_value:.4g}. This association does not establish predictive utility.")
    st.caption("All displayed results are generated from the SQLite database produced by the analysis pipeline.")


def sample_composition(sample_rows: pd.DataFrame) -> go.Figure:
    ordered = sample_rows.set_index("population").reindex(POPULATION_ORDER).reset_index()
    figure = go.Figure(go.Bar(
        y=ordered.population.map(POPULATION_LABELS), x=ordered.percentage,
        orientation="h", marker_color=CRIMSON, customdata=ordered["count"],
        text=ordered.percentage.map(lambda value: f"{value:.1f}%"), textposition="outside",
        hovertemplate="%{y}<br>%{x:.2f}% · %{customdata:,} cells<extra></extra>",
    ))
    chart_layout(figure, 290, dict(l=15, r=48, t=8, b=35))
    figure.update_yaxes(autorange="reversed", showgrid=False, title=None)
    figure.update_xaxes(title="Relative frequency (%)",
                        range=[0, max(40, ordered.percentage.max() * 1.2)],
                        showgrid=True, gridcolor=GRID)
    return figure


def population_distribution(frequencies: pd.DataFrame) -> go.Figure:
    figure = go.Figure()
    for population in POPULATION_ORDER:
        subset = frequencies.loc[frequencies.population == population]
        figure.add_trace(go.Box(
            name=POPULATION_LABELS[population], y=subset.percentage,
            line_color=CRIMSON, fillcolor=PALE_CRIMSON, marker_color=CRIMSON,
            boxpoints=False, showlegend=False,
            hovertemplate="%{y:.2f}%<extra>" + POPULATION_LABELS[population] + "</extra>",
        ))
    chart_layout(figure, 350, dict(l=45, r=15, t=12, b=55))
    figure.update_yaxes(title="Relative frequency (%)")
    return figure


def frequencies_page(database_revision: tuple[int, int]) -> None:
    page_header("02", "Immune Composition",
                "How is each sample distributed across the five measured immune-cell populations?")
    st.markdown("Relative frequency is the population count divided by the total count across "
                "the five measured populations within that sample.")
    frequencies = load_frequencies(database_revision)
    section_title("Sample composition", "Select a sample for its five-population profile")
    sample_id = st.selectbox("Sample ID", frequencies["sample"].drop_duplicates().tolist(),
                             index=0, help="Type in the selector to search by sample ID.")
    sample_rows = frequencies.loc[frequencies["sample"] == sample_id]
    total = int(sample_rows.total_count.iloc[0])
    st.markdown(f'<div class="sample-summary"><strong>{sample_id}</strong>'
                f'<span>{total:,} total measured cells</span></div>', unsafe_allow_html=True)
    st.plotly_chart(sample_composition(sample_rows), width="stretch", key="sample_composition")
    section_title("Across all samples", f"{frequencies['sample'].nunique():,} sample distributions")
    st.markdown("The spread of within-sample percentages shows how composition varies across the study.")
    st.plotly_chart(population_distribution(frequencies), width="stretch", key="population_distribution")
    section_title("Sample-level frequency data", "Complete Part 2 result")
    left, right = st.columns((1.2, 1), gap="medium")
    search = left.text_input("Find sample", placeholder="Search sample ID").strip()
    population = right.selectbox("Population", ("All populations", *POPULATION_ORDER),
                                  format_func=lambda value: POPULATION_LABELS.get(value, value))
    filtered = frequencies
    if population != "All populations":
        filtered = filtered.loc[filtered.population == population]
    if search:
        filtered = filtered.loc[filtered["sample"].str.contains(search, case=False, regex=False)]
    page_size = 250
    page_count = max(1, math.ceil(len(filtered) / page_size))
    st.caption(f"{len(filtered):,} rows · {filtered['sample'].nunique():,} samples")
    page = st.number_input("Page", min_value=1, max_value=page_count, value=1, step=1)
    page_rows = filtered.iloc[(page - 1) * page_size:page * page_size].copy()
    page_rows["population"] = page_rows["population"].map(POPULATION_LABELS)
    st.dataframe(page_rows, hide_index=True, width="stretch", height=360,
                 column_config={
                     "sample": st.column_config.TextColumn("Sample"),
                     "total_count": st.column_config.NumberColumn("Total count", format="%d"),
                     "population": st.column_config.TextColumn("Population"),
                     "count": st.column_config.NumberColumn("Count", format="%d"),
                     "percentage": st.column_config.NumberColumn("Percentage", format="%.2f%%"),
                 })
    st.download_button("Download complete frequency table", frequencies.to_csv(index=False),
                       file_name="frequency_table.csv", mime="text/csv")


def response_boxplots(cohort: pd.DataFrame, results: pd.DataFrame) -> go.Figure:
    figure = make_subplots(rows=2, cols=3, shared_yaxes=True, vertical_spacing=.22,
                           horizontal_spacing=.09,
                           subplot_titles=[POPULATION_LABELS[population] for population in POPULATION_ORDER])
    results_by_population = results.set_index("population")
    for position, population in enumerate(POPULATION_ORDER):
        row, column = divmod(position, 3)
        subset = cohort.loc[cohort.population == population]
        for response, label in (("no", "Nonresponders"), ("yes", "Responders")):
            group = subset.loc[subset.response == response]
            figure.add_trace(go.Box(
                y=group.percentage, name=label, legendgroup=response,
                showlegend=position == 0, line_color=GROUP_COLORS[response],
                fillcolor=GROUP_COLORS[response], opacity=.72, boxpoints=False,
                marker_color=GROUP_COLORS[response],
                hovertemplate="%{y:.2f}%<extra>" + label + "</extra>",
            ), row=row + 1, col=column + 1)
        statistic = results_by_population.loc[population]
        axis_name = "" if position == 0 else str(position + 1)
        figure.add_annotation(xref=f"x{axis_name} domain", yref=f"y{axis_name} domain",
                              x=.5, y=.92, showarrow=False,
                              text=f"p = {statistic.raw_p_value:.3g} · q = {statistic.fdr_q_value:.3g}",
                              font=dict(size=10, color=MUTED))
    chart_layout(figure, 670, dict(l=48, r=18, t=92, b=35))
    figure.update_layout(boxmode="group", font=dict(color=CHARCOAL, size=11))
    figure.update_yaxes(range=[0, max(50, cohort.percentage.max() * 1.05)])
    figure.update_yaxes(title="Relative frequency (%)", row=1, col=1)
    figure.update_yaxes(title="Relative frequency (%)", row=2, col=1)
    for annotation in figure.layout.annotations[:5]:
        annotation.update(font=dict(size=13, color=CHARCOAL))
    return figure


def display_results(table: pd.DataFrame) -> None:
    display = table[["population", "responder_median", "nonresponder_median",
                     "effect_size", "raw_p_value", "fdr_q_value",
                     "nominal_significant", "fdr_significant"]].copy()
    display["population"] = display["population"].map(POPULATION_LABELS)
    display["interpretation"] = display.apply(
        lambda row: "FDR significant" if row.fdr_significant else
                    "Nominal only" if row.nominal_significant else "No nominal evidence", axis=1)
    display = display.drop(columns=["nominal_significant", "fdr_significant"])
    def shade(row):
        color = "background-color: #EBD4DA" if row["fdr_q_value"] < .05 else (
            "background-color: #F7ECEF" if row["raw_p_value"] < .05 else "")
        return [color] * len(row)
    st.dataframe(
        display.style.apply(shade, axis=1), hide_index=True, width="stretch",
        column_config={
            "population": st.column_config.TextColumn("Population"),
            "responder_median": st.column_config.NumberColumn("Responder median (%)", format="%.2f"),
            "nonresponder_median": st.column_config.NumberColumn("Nonresponder median (%)", format="%.2f"),
            "effect_size": st.column_config.NumberColumn("Rank-biserial", format="%.3f"),
            "raw_p_value": st.column_config.NumberColumn("Raw p", format="%.4f"),
            "fdr_q_value": st.column_config.NumberColumn("FDR q", format="%.4f"),
            "interpretation": st.column_config.TextColumn("Interpretation"),
        },
    )


def response_page(database_revision: tuple[int, int]) -> None:
    page_header("03", "Treatment Response",
                "Do immune-cell frequencies differ between patients who respond to miraclib and those who do not?")
    analysis = load_response_analysis(database_revision)
    summary = analysis["summary"]
    st.markdown('<div class="cohort-line">Melanoma · miraclib · PBMC</div>', unsafe_allow_html=True)
    metrics([("Pooled sample observations", f"{summary['samples']:,}"),
             ("Unique subjects", f"{summary['subjects']:,}"),
             ("Responder samples", f"{summary['responder_samples']:,}"),
             ("Nonresponder samples", f"{summary['nonresponder_samples']:,}")])
    section_title("Key finding")
    lead = strongest(analysis["pooled"])
    st.markdown(f'**Strongest observed difference: {POPULATION_LABELS[lead.population]}**')
    metrics([("Responder median", f"{lead.responder_median:.2f}%"),
             ("Nonresponder median", f"{lead.nonresponder_median:.2f}%"),
             ("Rank-biserial", f"{lead.effect_size:.3f}"),
             ("Raw p", f"{lead.raw_p_value:.3g}"),
             ("FDR q", f"{lead.fdr_q_value:.3g}")])
    st.markdown(evidence_text(lead))
    section_title("Five-population comparison", "Pooled samples · two-sided Mann–Whitney U")
    st.caption("Boxes show response-group distributions; p is unadjusted and q controls FDR across five populations.")
    st.plotly_chart(response_boxplots(analysis["cohort"], analysis["pooled"]),
                    width="stretch", key="response_boxplots")
    section_title("Statistical evidence")
    display_results(analysis["pooled"])
    st.caption("Responder and nonresponder distributions were compared using two-sided Mann–Whitney U tests. "
               "Benjamini–Hochberg FDR adjustment covers the five populations. "
               "Pooled longitudinal samples from one subject are not fully independent.")
    section_title("Robustness check", "Subject-level and timepoint sensitivity analyses")
    st.markdown("**Why check this?** The primary comparison pools three samples from each patient, "
                "so the same person is counted three times. These checks count each patient once: "
                "as a mean across timepoints, as a single Day 0 sample, or as a change from Day 0 to Day 14.")
    subject_lead = strongest(analysis["subject_mean"])
    st.markdown(f"In the subject-mean comparison, **{POPULATION_LABELS[subject_lead.population]}** "
                f"had the smallest raw p (p = {subject_lead.raw_p_value:.3g}; "
                f"FDR q = {subject_lead.fdr_q_value:.3g}). {evidence_text(subject_lead)}")
    change = analysis["subject_change"]
    significant_change = change.loc[change.fdr_significant]
    if significant_change.empty:
        st.markdown("No population's change from Day 0 to Day 14 met the FDR threshold (q < 0.05).")
    else:
        for row in significant_change.itertuples():
            st.markdown(f"**Change from Day 0 to Day 14: {POPULATION_LABELS[row.population]}**")
            metrics([("Responder median change", f"{row.responder_median:+.2f} pp"),
                     ("Nonresponder median change", f"{row.nonresponder_median:+.2f} pp"),
                     ("Raw p", f"{row.raw_p_value:.3g}"),
                     ("FDR q", f"{row.fdr_q_value:.3g}")])
            st.markdown(f"{POPULATION_LABELS[row.population]} changed differently in responders and "
                        "nonresponders during treatment, and this difference passed FDR correction.")
        st.caption("Treat this as exploratory. FDR correction covers the five populations within each "
                   "analysis, not all four analyses together, and a change during treatment cannot "
                   "predict response before treatment starts.")
    with st.expander("Subject-mean results"):
        display_results(analysis["subject_mean"])
    with st.expander("Day 0 results"):
        st.caption("One day-zero observation per subject. Day 0 is treated as baseline, "
                   "but pretreatment collection is not confirmed by the source data.")
        display_results(analysis["baseline"])
    with st.expander("Day 14 − day 0 change results"):
        st.caption("Each subject contributes one percentage-point change per population. "
                   "This is a longitudinal association, not a baseline prediction.")
        display_results(analysis["subject_change"])


def distribution_bars(rows: list[dict], label_key: str, count_key: str,
                      labels: dict[str, str], color: str) -> None:
    total = sum(row[count_key] for row in rows)
    for row in rows:
        label = labels.get(row[label_key], row[label_key])
        count = row[count_key]
        share = count / total * 100 if total else 0
        st.markdown(f'<div class="dist-label"><span>{label}</span>'
                    f'<strong>{count:,} <small>({share:.1f}%)</small></strong></div>'
                    f'<div class="dist-track"><div style="width:{share:.1f}%;background:{color}"></div></div>',
                    unsafe_allow_html=True)


def baseline_page(database_revision: tuple[int, int]) -> None:
    page_header("04", "Cohort Explorer", "Who makes up the baseline melanoma cohort?")
    result = load_baseline_subset(database_revision)
    cohort = result["cohort"]
    tables = result["tables"]
    st.markdown('<div class="cohort-line">Melanoma · miraclib · PBMC · time = 0</div>', unsafe_allow_html=True)
    metrics([("Qualifying samples", f"{len(cohort):,}"),
             ("Unique subjects", f"{len({row['subject'] for row in cohort}):,}")])
    section_title("Baseline cohort samples", f"{len(cohort):,} samples")
    st.markdown("All melanoma PBMC samples at time 0 from patients treated with miraclib.")
    cohort_table = pd.DataFrame(cohort)[["sample", "subject", "project", "response", "sex"]]
    st.dataframe(cohort_table, hide_index=True, width="stretch", height=320,
                 column_config={
                     "sample": st.column_config.TextColumn("Sample"),
                     "subject": st.column_config.TextColumn("Subject"),
                     "project": st.column_config.TextColumn("Project"),
                     "response": st.column_config.TextColumn("Response"),
                     "sex": st.column_config.TextColumn("Sex"),
                 })
    st.download_button("Download baseline cohort", cohort_table.to_csv(index=False),
                       file_name="part4_cohort.csv", mime="text/csv")
    section_title("Cohort composition", "Counts and shares at the correct observational level")
    st.markdown("**Samples by project** · unique samples")
    distribution_bars(tables["samples_by_project"], "project", "sample_count", {}, CRIMSON)
    st.markdown("**Subjects by response** · unique subjects")
    distribution_bars(tables["subjects_by_response"], "response", "subject_count",
                      {"yes": "Responders", "no": "Nonresponders"}, CRIMSON)
    st.markdown("**Subjects by sex** · unique subjects")
    distribution_bars(tables["subjects_by_sex"], "sex", "subject_count",
                      {"M": "Male", "F": "Female"}, DARK_CRIMSON)
    section_title("Additional subset query")
    st.markdown("**Among male melanoma responders at time 0, across all sample and treatment types, "
                "what is the average raw B-cell count?**")
    query = result["b_cell_query"]
    metrics([("Average raw B-cell count", f"{query['average_b_cell_count']:,.2f} cells"),
             ("Qualifying samples", f"{query['qualifying_samples']:,}"),
             ("Unique subjects", f"{query['qualifying_subjects']:,}")])
    st.caption(f"Sample types represented: {', '.join(query['sample_types'])} · "
               f"Treatments represented: {', '.join(query['treatments'])}. "
               "No sample-type or treatment filter is applied to this query.")


def main() -> None:
    st.set_page_config(page_title="Immune Cell Analysis", layout="wide",
                       initial_sidebar_state="expanded")
    st.markdown("""
    <style>
      :root {color-scheme: light;}
      .stApp {background: #F7F5F3; color: #202124;}
      .block-container {max-width: 1160px; padding-top: 5rem; padding-bottom: 3rem;}
      h1 {font-size: 2.25rem !important; letter-spacing: -.035em; line-height: 1.16; margin: .15rem 0 .55rem !important;}
      h2 {font-size: 1.13rem !important; letter-spacing: -.015em;}
      p {line-height: 1.55;}
      section[data-testid="stSidebar"] {background: #FFFFFF; border-right: 1px solid #DED9D6;}
      section[data-testid="stSidebar"] > div {padding-top: 2rem;}
      .sidebar-brand {font-size: .86rem; font-weight: 800; letter-spacing: .095em; color: #65192A;}
      .sidebar-sub {font-size: .78rem; color: #6E6A69; margin: .15rem 0 1.7rem;}
      div[data-testid="stRadio"] > label {display:none;}
      div[data-testid="stRadio"] div[role="radiogroup"] {gap: .18rem;}
      div[data-testid="stRadio"] label {border-left: 3px solid transparent; padding: .55rem .65rem; border-radius: 0; width: 100%;}
      div[data-testid="stRadio"] label:has(input:checked) {border-left-color: #8F263D; background: #F7ECEF; color: #65192A; font-weight: 650;}
      label[data-testid="stRadioOption"] > div > div:first-child {display:none !important;}
      .eyebrow {font-size: .73rem; font-weight: 800; letter-spacing: .16em; color: #8F263D;}
      .page-question {font-size: 1rem; color: #6E6A69; margin: 0 0 1.2rem;}
      .page-rule {height: 2px; background: #DED9D6; border-left: 62px solid #8F263D; margin-bottom: 1.6rem;}
      .section-heading {display:flex; align-items:baseline; justify-content:space-between; gap:1rem; border-top:1px solid #DED9D6; margin: 2.15rem 0 .9rem; padding-top: .75rem;}
      .section-heading h2 {margin:0; font-size:1.18rem !important; font-weight:700 !important;}
      .section-heading span {font-size:.75rem; color:#6E6A69; text-align:right;}
      .report-metric {border-top:2px solid #8F263D; padding:.55rem 0 .45rem; min-height:4.5rem;}
      .metric-value {font-size:1.45rem; font-weight:700; color:#202124; letter-spacing:-.025em; white-space:nowrap;}
      .metric-label {font-size:.76rem; color:#6E6A69; margin-top:.12rem;}
      .structure-step {border-top:1px solid #DED9D6; padding:.8rem 0; display:flex; flex-direction:column;}
      .structure-step strong {font-size:1.35rem; color:#65192A;}
      .structure-step span {font-size:.78rem; color:#6E6A69;}
      .question-row {display:flex; align-items:baseline; gap:1.1rem; border-top:1px solid #DED9D6; padding:.56rem 0;}
      .question-row span {font-size:.75rem; font-weight:800; color:#8F263D;}
      .question-row p {margin:0;}
      .cohort-line {font-size:.85rem; color:#65192A; font-weight:700; letter-spacing:.02em; margin-bottom:1rem;}
      .sample-summary {display:flex; justify-content:space-between; border-bottom:1px solid #DED9D6; padding:.15rem 0 .6rem;}
      .sample-summary span {color:#6E6A69;}
      .dist-label {display:flex; justify-content:space-between; margin:.55rem 0 .25rem; font-size:.85rem;}
      .dist-label small {font-weight:400; color:#6E6A69;}
      .dist-track {height:.5rem; background:#EAE6E3; margin-bottom:.8rem;}
      .dist-track div {height:100%;}
      div[data-testid="stDataFrame"] {border:1px solid #DED9D6;}
      a {color:#8F263D;}
      @media (max-width: 760px) {
        .block-container {padding-top:4.5rem;}
        .metric-value {font-size:1.13rem;}
        .section-heading {display:block;}
        .section-heading span {display:block; text-align:left; margin-top:.2rem;}
      }
    </style>
    """, unsafe_allow_html=True)
    st.sidebar.markdown('<div class="sidebar-brand">IMMUNE CELL ANALYSIS</div>'
                        '<div class="sidebar-sub">Technical assessment</div>',
                        unsafe_allow_html=True)
    section = st.sidebar.radio("Sections", PAGES, label_visibility="collapsed")
    if not DB_PATH.is_file():
        st.error("Database not found. Run `make pipeline` before starting the dashboard.")
        st.stop()
    stat = DB_PATH.stat()
    database_revision = (stat.st_mtime_ns, stat.st_size)
    try:
        if section == PAGES[0]:
            overview_page(database_revision)
        elif section == PAGES[1]:
            frequencies_page(database_revision)
        elif section == PAGES[2]:
            response_page(database_revision)
        else:
            baseline_page(database_revision)
    except (sqlite3.Error, ValueError) as exc:
        st.error(f"Could not load the analysis: {exc}")
        st.stop()


if __name__ == "__main__":
    main()