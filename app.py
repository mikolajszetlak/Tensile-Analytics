import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import linregress
import plotly.graph_objects as go
import io

# Konfiguracja strony
st.set_page_config(
    page_title="Tensile Analytics Pro | ISO 6892-1 / ASTM E8",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- SILNIK OBLICZENIOWY ---
def process_tensile_dataset(
    df_raw: pd.DataFrame,
    force_col: str,
    disp_col: str,
    force_unit: str,
    disp_unit: str,
    l0: float,
    s0: float,
    sheet_name: str
):
    """
    Deterministyczny silnik wyznaczania parametrów wytrzymałościowych (ISO 6892-1 / ASTM E8M).
    Zawiera filtrację uślizgu uchwytów, dynamiczny fallback i zabezpieczenie numeryczne.
    """
    # 1. Konwersja numeryczna i czyszczenie
    force_series = pd.to_numeric(df_raw[force_col], errors='coerce')
    disp_series = pd.to_numeric(df_raw[disp_col], errors='coerce')
    
    clean_mask = force_series.notna() & disp_series.notna()
    force_arr = force_series[clean_mask].to_numpy(dtype=np.float64)
    disp_arr = disp_series[clean_mask].to_numpy(dtype=np.float64)
    
    if len(force_arr) < 20:
        return None

    # 2. Standaryzacja do układu bazowego: Siła [N], Przemieszczenie [mm]
    if force_unit == "kN":
        force_arr = force_arr * 1e3
    elif force_unit == "MN":
        force_arr = force_arr * 1e6
        
    if disp_unit in ["μm", "um"]:
        disp_arr = disp_arr / 1e3
    elif disp_unit == "m":
        disp_arr = disp_arr * 1e3

    # 3. Wyznaczenie wektorów naprężeń sigma [MPa] i odkształceń epsilon [-]
    stress = force_arr / s0
    strain = disp_arr / l0
    
    # 4. Parametry wytrzymałościowe graniczne (Rm, Ru, A)
    rm_idx = int(np.argmax(stress))
    if rm_idx < 5:
        rm_idx = len(stress) - 1
        
    fm_kn = float(force_arr[rm_idx] / 1000.0)
    rm_mpa = float(stress[rm_idx])
    
    fu_kn = float(force_arr[-1] / 1000.0)
    ru_mpa = float(stress[-1])
    a_pct = float(strain[-1] * 100.0)
    
    # 5. Ograniczenie przestrzeni poszukiwań strefy Hooke'a (odcięcie uślizgu uchwytów)
    stress_min_cutoff = 0.15 * rm_mpa
    stress_max_cutoff = 0.75 * rm_mpa
    
    search_indices = np.where(
        (stress >= stress_min_cutoff) & 
        (stress <= stress_max_cutoff) & 
        (np.arange(len(stress)) <= rm_idx)
    )[0]
    
    if len(search_indices) < 15:
        search_indices = np.arange(0, max(15, int(rm_idx * 0.75)))

    # 6. Adaptacyjne ruchome okno szukające bezwzględnego maksimum sztywności
    window_size = max(6, min(50, int(len(search_indices) * 0.25)))
    best_slope = 0.0
    best_intercept = 0.0
    best_r2 = 0.0
    best_idx = search_indices[0]
    
    for start_pos in range(0, len(search_indices) - window_size + 1):
        cur_idx_window = search_indices[start_pos : start_pos + window_size]
        w_strain = strain[cur_idx_window]
        w_stress = stress[cur_idx_window]
        
        if np.ptp(w_strain) == 0:
            continue
            
        slope, intercept, r_val, _, _ = linregress(w_strain, w_stress)
        r2 = r_val ** 2
        
        if r2 > 0.985 and slope > best_slope:
            best_slope = slope
            best_intercept = intercept
            best_r2 = r2
            best_idx = cur_idx_window[0]

    # Fallback dla danych zaszumionych
    if best_slope <= 0.0:
        for start_pos in range(0, len(search_indices) - window_size + 1):
            cur_idx_window = search_indices[start_pos : start_pos + window_size]
            w_strain = strain[cur_idx_window]
            w_stress = stress[cur_idx_window]
            if np.ptp(w_strain) == 0:
                continue
            slope, intercept, r_val, _, _ = linregress(w_strain, w_stress)
            r2 = r_val ** 2
            if slope > 0 and r2 > best_r2:
                best_slope = slope
                best_intercept = intercept
                best_r2 = r2
                best_idx = cur_idx_window[0]

    # Bezpiecznik przed dzieleniem przez zero
    if best_slope <= 0.0:
        delta_e = strain[rm_idx] - strain[0]
        best_slope = (stress[rm_idx] - stress[0]) / delta_e if delta_e > 0 else 1000.0
        best_intercept = stress[0] - best_slope * strain[0]

    # 7. Kompensacja luzu początkowego (Toe Compensation)
    toe_offset = -best_intercept / best_slope
    strain_corr = strain - toe_offset
    
    # 8. Wyznaczenie Rp0.2 (Offset 0.2%)
    sigma_offset = best_slope * (strain_corr - 0.002)
    diff = stress - sigma_offset
    
    min_eval_idx = max(best_idx + window_size, 0)
    valid_range = np.where((strain_corr >= 0.002) & (np.arange(len(strain)) >= min_eval_idx))[0]
    
    rp02_val = None
    rp02_strain = None
    
    if len(valid_range) > 1:
        for idx in valid_range[:-1]:
            if diff[idx] * diff[idx + 1] <= 0:
                denom = abs(diff[idx]) + abs(diff[idx + 1])
                weight = abs(diff[idx]) / denom if denom > 0 else 0.5
                rp02_strain = float(strain_corr[idx] + weight * (strain_corr[idx + 1] - strain_corr[idx]))
                rp02_val = float(stress[idx] + weight * (stress[idx + 1] - stress[idx]))
                break
                
    return {
        "sheet": str(sheet_name),
        "Fm_kN": fm_kn,
        "Rm_MPa": rm_mpa,
        "Rm_strain_pct": float(strain_corr[rm_idx] * 100.0),
        "Rp02_MPa": rp02_val,
        "Rp02_strain_pct": (rp02_strain * 100.0) if rp02_strain else None,
        "Fu_kN": fu_kn,
        "Ru_MPa": ru_mpa,
        "A_pct": a_pct,
        "strain_corr_pct": strain_corr * 100.0,
        "stress": stress,
        "slope": best_slope,
        "toe_offset_pct": toe_offset * 100.0
    }

# --- INTERFEJS UŻYTKOWNIKA ---
st.title("🔬 Moduł Zaawansowanej Analizy Próby Rozciągania")
st.caption("Standardy metrologiczne: PN-EN ISO 6892-1 / ASTM E8M (Deterministyczna Kompensacja Luzu)")

with st.sidebar:
    st.header("⚙️ 1. Plik Źródłowy")
    uploaded_file = st.file_uploader("Załaduj arkusz pomiarowy (.xlsx, .xls)", type=["xlsx", "xls"])

if not uploaded_file:
    st.info("👆 Załaduj plik Excel z wynikami próby rozciągania w lewym panelu, aby rozpocząć pracę.")
else:
    xl_doc = pd.ExcelFile(uploaded_file)
    available_sheets = xl_doc.sheet_names
    
    # SEKCJA INSPEKCJI SUROWYCH WIERSZY
    with st.expander("🔍 Podgląd Surowego Nagłówka i Metadanych Próbki (Wskaźnik Wierszy)", expanded=True):
        col_insp1, col_insp2 = st.columns([3, 1])
        with col_insp1:
            preview_sheet = st.selectbox("Wybierz arkusz do inspekcji nagłówka:", available_sheets, key="insp_sheet")
        with col_insp2:
            n_rows_preview = st.number_input("Liczba wierszy do wyświetlenia:", min_value=3, max_value=30, value=10, step=1)
            
        # Wczytanie surowych wierszy bez nagłówka
        df_raw_preview = pd.read_excel(uploaded_file, sheet_name=preview_sheet, header=None, nrows=int(n_rows_preview))
        
        # Formatowanie indeksu jako "Wiersz X" dla jednoznacznej identyfikacji
        df_raw_preview.index = [f"Wiersz {i}" for i in range(len(df_raw_preview))]
        st.dataframe(df_raw_preview, use_container_width=True)
        st.caption("💡 Sprawdź powyżej: w którym wierszu znajdują się $L_0$ i $S_0$ oraz od którego wiersza zaczynają się nazwy kolumn.")

    with st.sidebar:
        st.markdown("---")
        st.header("📐 2. Geometria Próbki")
        col_g1, col_g2 = st.columns(2)
        with col_g1:
            manual_l0 = st.number_input("Baza L₀ [mm]", value=40.0, step=1.0, min_value=1.0)
        with col_g2:
            manual_s0 = st.number_input("Przekrój S₀ [mm²]", value=73.125, step=0.5, min_value=0.1)

        st.markdown("---")
        st.header("🛠️ 3. Konfiguracja Struktury Danych")
        header_row_idx = st.number_input(
            "Indeks wiersza z nagłówkami kolumn (np. 5):",
            value=5,
            min_value=0,
            max_value=50,
            step=1,
            help="Numer wiersza z nazwami kolumn odczytany z tabeli inspekcji powyżej."
        )
        
        skip_after_header = st.number_input(
            "Pomiń wiersze pod nagłówkiem (np. jednostki):",
            value=1,
            min_value=0,
            max_value=10,
            step=1,
            help="Liczba wierszy z jednostkami bezpośrednio pod nagłówkiem (np. s, N, mm)."
        )

    # Odczyt roboczy arkusza po konfiguracji nagłówków
    df_preview = pd.read_excel(uploaded_file, sheet_name=available_sheets[0], header=int(header_row_idx))
    if skip_after_header > 0:
        df_preview = df_preview.iloc[int(skip_after_header):].reset_index(drop=True)
        
    df_preview.columns = [str(col).strip() for col in df_preview.columns]
    column_options = df_preview.columns.tolist()

    with st.sidebar:
        st.markdown("---")
        st.header("🎯 4. Mapowanie Kolumn i Jednostek")
        default_force_idx = 1 if len(column_options) > 1 else 0
        default_disp_idx = 2 if len(column_options) > 2 else 0
        
        force_column = st.selectbox("Kolumna Siły (F):", column_options, index=default_force_idx)
        force_unit = st.selectbox("Jednostka Siły w pliku:", ["N", "kN", "MN"], index=0)
        
        disp_column = st.selectbox("Kolumna Przemieszczenia (ΔL):", column_options, index=default_disp_idx)
        disp_unit = st.selectbox("Jednostka Przemieszczenia w pliku:", ["mm", "μm", "m"], index=0)

        st.markdown("---")
        st.header("📊 5. Tryb Analizy")
        analysis_mode = st.radio(
            "Wybierz tryb generowania raportu:",
            ["Pojedyncza Próbka (Raport Szczegółowy)", "Zbiorcza Seria Pomiarowa (Multi-Sheet)"]
        )

    # =========================================================================
    # TRYB 1: POJEDYNCZA PRÓBKA (RAPORT SZCZEGÓŁOWY)
    # =========================================================================
    if analysis_mode == "Pojedyncza Próbka (Raport Szczegółowy)":
        with st.sidebar:
            selected_sheet = st.selectbox("Wybierz arkusz próbki do analizy:", available_sheets)
            
        df_target = pd.read_excel(uploaded_file, sheet_name=selected_sheet, header=int(header_row_idx))
        if skip_after_header > 0:
            df_target = df_target.iloc[int(skip_after_header):].reset_index(drop=True)
        df_target.columns = [str(col).strip() for col in df_target.columns]
        
        res = process_tensile_dataset(
            df_target, force_column, disp_column, force_unit, disp_unit, manual_l0, manual_s0, selected_sheet
        )
        
        if res is None:
            st.error("Błąd: Wybrany arkusz nie zawiera wystarczającej liczby poprawnych punktów pomiarowych.")
        else:
            col1, col2, col3, col4, col5 = st.columns(5)
            col1.metric("Wytrzymałość Doraźna Rm", f"{res['Rm_MPa']:.1f} MPa")
            col2.metric("Granica Plastyczności Rp0.2", f"{res['Rp02_MPa']:.1f} MPa" if res['Rp02_MPa'] else "Brak")
            col3.metric("Naprężenie Zerwania Ru", f"{res['Ru_MPa']:.1f} MPa")
            col4.metric("Wydłużenie Całkowite A", f"{res['A_pct']:.1f} %")
            col5.metric("Maksymalna Siła Fm", f"{res['Fm_kN']:.2f} kN")
            
            fig = go.Figure()
            
            # Krzywa rozciągania
            fig.add_trace(go.Scatter(
                x=res["strain_corr_pct"],
                y=res["stress"],
                mode='lines',
                name=f'Krzywa σ-ε (Próbka {selected_sheet})',
                line=dict(color='#0066cc', width=2.5)
            ))
            
            # Ograniczona prosta offsetowa Rp0.2
            rp02_target = res["Rp02_MPa"] if res["Rp02_MPa"] else (res["Rm_MPa"] * 0.8)
            max_offset_stress = min(res["Rm_MPa"] * 0.95, rp02_target * 1.15)
            max_offset_strain = (max_offset_stress / res["slope"])
            
            e_range = np.linspace(0, max_offset_strain, 50)
            fig.add_trace(go.Scatter(
                x=(e_range + 0.002) * 100,
                y=res["slope"] * e_range,
                mode='lines',
                name='Offset 0.2%',
                line=dict(color='rgba(230, 0, 0, 0.7)', dash='dash', width=1.8)
            ))
            
            # Punkt Rp0.2
            if res["Rp02_MPa"] and res["Rp02_strain_pct"]:
                fig.add_trace(go.Scatter(
                    x=[res["Rp02_strain_pct"]],
                    y=[res["Rp02_MPa"]],
                    mode='markers+text',
                    name=f'Rp0.2 = {res["Rp02_MPa"]:.1f} MPa',
                    text=[f"<b>Rp0.2:</b> {res['Rp02_MPa']:.1f} MPa"],
                    textposition="top left",
                    marker=dict(color='#d90429', size=9, symbol='diamond')
                ))
                
            # Punkt Rm
            fig.add_trace(go.Scatter(
                x=[res["Rm_strain_pct"]],
                y=[res["Rm_MPa"]],
                mode='markers+text',
                name=f'Rm = {res["Rm_MPa"]:.1f} MPa',
                text=[f"<b>Rm:</b> {res['Rm_MPa']:.1f} MPa"],
                textposition="top center",
                marker=dict(color='#f77f00', size=10, symbol='circle')
            ))
            
            # Punkt Ru
            fig.add_trace(go.Scatter(
                x=[res["strain_corr_pct"][-1]],
                y=[res["Ru_MPa"]],
                mode='markers+text',
                name=f'Ru = {res["Ru_MPa"]:.1f} MPa',
                text=[f"<b>Ru:</b> {res['Ru_MPa']:.1f} MPa (A = {res['A_pct']:.1f}%)"],
                textposition="bottom left",
                marker=dict(color='#2b2d42', size=9, symbol='x')
            ))
            
            y_max_limit = float(res["Rm_MPa"] * 1.12)
            x_max_limit = float(res["strain_corr_pct"][-1] * 1.08)
            
            fig.update_layout(
                title=dict(
                    text=f"<b>Charakterystyka Rozciągania – Próbka {selected_sheet}</b>", 
                    font=dict(size=17, color="#f8f9fa")
                ),
                xaxis=dict(
                    title="Odkształcenie skorygowane ε_corr [%]",
                    range=[-0.5, x_max_limit],
                    gridcolor="rgba(255, 255, 255, 0.1)",
                    zeroline=True,
                    zerolinecolor="rgba(255, 255, 255, 0.2)"
                ),
                yaxis=dict(
                    title="Naprężenie inżynierskie σ [MPa]",
                    range=[0, y_max_limit],
                    gridcolor="rgba(255, 255, 255, 0.1)",
                    zeroline=True,
                    zerolinecolor="rgba(255, 255, 255, 0.2)"
                ),
                hovermode="closest",
                template="plotly_dark",
                legend=dict(
                    yanchor="top", 
                    y=0.95, 
                    xanchor="left", 
                    x=0.03, 
                    bgcolor="rgba(20, 20, 20, 0.75)",
                    bordercolor="rgba(255, 255, 255, 0.15)",
                    borderwidth=1
                ),
                height=580,
                margin=dict(l=60, r=40, t=60, b=50)
            )
            
            st.plotly_chart(fig, use_container_width=True)

    # =========================================================================
    # TRYB 2: ZBIORCZA SERIA POMIAROWA (MULTI-SHEET)
    # =========================================================================
    else:
        with st.sidebar:
            selected_sheets = st.multiselect(
                "Wybierz arkusze do porównania zbiorczego:",
                available_sheets,
                default=available_sheets
            )
            
        if not selected_sheets:
            st.warning("Zaznacz co najmniej jeden arkusz do wygenerowania raportu serii.")
        else:
            batch_results = []
            fig_cum = go.Figure()
            max_rm_overall = 0.0
            
            for sheet in selected_sheets:
                df_item = pd.read_excel(uploaded_file, sheet_name=sheet, header=int(header_row_idx))
                if skip_after_header > 0:
                    df_item = df_item.iloc[int(skip_after_header):].reset_index(drop=True)
                df_item.columns = [str(col).strip() for col in df_item.columns]
                
                res_item = process_tensile_dataset(
                    df_item, force_column, disp_column, force_unit, disp_unit, manual_l0, manual_s0, sheet
                )
                
                if res_item is not None:
                    if res_item["Rm_MPa"] > max_rm_overall:
                        max_rm_overall = res_item["Rm_MPa"]
                        
                    batch_results.append({
                        "Próbka": f"Próbka {sheet}",
                        "Fm [kN]": round(res_item["Fm_kN"], 2),
                        "Rm [MPa]": round(res_item["Rm_MPa"], 1),
                        "Rp0.2 [MPa]": round(res_item["Rp02_MPa"], 1) if res_item["Rp02_MPa"] else np.nan,
                        "Fu [kN]": round(res_item["Fu_kN"], 2),
                        "Ru [MPa]": round(res_item["Ru_MPa"], 1),
                        "A [%]": round(res_item["A_pct"], 1)
                    })
                    
                    fig_cum.add_trace(go.Scatter(
                        x=res_item["strain_corr_pct"],
                        y=res_item["stress"],
                        mode='lines',
                        name=f'Próbka {sheet}',
                        line=dict(width=1.8)
                    ))
                    
            if not batch_results:
                st.error("Żaden z wybranych arkuszy nie zawierał poprawnych danych.")
            else:
                df_batch = pd.DataFrame(batch_results)
                
                st.subheader("📈 Skumulowany Przebieg Krzywych Rozciągania Serii")
                fig_cum.update_layout(
                    xaxis=dict(
                        title="Odkształcenie skorygowane ε_corr [%]",
                        gridcolor="rgba(255, 255, 255, 0.1)"
                    ),
                    yaxis=dict(
                        title="Naprężenie inżynierskie σ [MPa]",
                        range=[0, max_rm_overall * 1.12],
                        gridcolor="rgba(255, 255, 255, 0.1)"
                    ),
                    template="plotly_dark",
                    height=540,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                )
                st.plotly_chart(fig_cum, use_container_width=True)
                
                st.subheader("📋 Zestawienie Wyników Pomiarowych Próbek")
                st.dataframe(df_batch, use_container_width=True)
                
                st.subheader("📊 Analiza Statystyczna Serii Badawczej")
                stat_cols = ["Rm [MPa]", "Rp0.2 [MPa]", "Ru [MPa]", "A [%]"]
                
                stats_data = []
                for col in stat_cols:
                    mean_v = df_batch[col].mean()
                    std_v = df_batch[col].std()
                    cov_v = (std_v / mean_v * 100.0) if mean_v != 0 else 0.0
                    stats_data.append({
                        "Parametr": col,
                        "Średnia (μ)": f"{mean_v:.2f}",
                        "Odchylenie std (σ)": f"{std_v:.2f}",
                        "Współczynnik zmienności V [%]": f"{cov_v:.2f}%"
                    })
                    
                st.table(pd.DataFrame(stats_data))
                
                csv_buffer = io.StringIO()
                df_batch.to_csv(csv_buffer, index=False, sep=';')
                st.download_button(
                    label="📥 Pobierz Pełny Raport CSV",
                    data=csv_buffer.getvalue(),
                    file_name="raport_serii_rozciagania.csv",
                    mime="text/csv"
                )