import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import linregress
import plotly.graph_objects as go
import io
import re

# Konfiguracja strony
st.set_page_config(
    page_title="Tensile Analytics  | SZETLAK",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- UNIWERSALNA INGESTIA PLIKÓW (EXCEL & CSV) ---
def detect_csv_separator(file_obj) -> str:
    """Wykrywa separator kolumn w pliku CSV na podstawie analizy pierwszych wierszy."""
    file_obj.seek(0)
    sample_bytes = file_obj.read(4096)
    file_obj.seek(0)
    
    try:
        sample_text = sample_bytes.decode('utf-8')
    except UnicodeDecodeError:
        sample_text = sample_bytes.decode('cp1250', errors='ignore')
        
    lines = [line for line in sample_text.splitlines() if line.strip()][:10]
    if not lines:
        return ';'
        
    counts = {
        ';': sum(line.count(';') for line in lines),
        ',': sum(line.count(',') for line in lines),
        '\t': sum(line.count('\t') for line in lines)
    }
    detected = max(counts, key=counts.get)
    return detected if counts[detected] > 0 else ';'

def load_raw_dataframe(file_obj, item_name, header_row=None, nrows=None, skip_after=0, manual_sep=None):
    """Uniwersalny parser odczytujący wycinek danych z pliku Excel lub CSV."""
    file_name = file_obj.name.lower()
    
    if file_name.endswith(('.xlsx', '.xls')):
        df = pd.read_excel(file_obj, sheet_name=item_name, header=header_row, nrows=nrows)
        if header_row is not None and skip_after > 0:
            df = df.iloc[int(skip_after):].reset_index(drop=True)
        return df
        
    elif file_name.endswith('.csv'):
        file_obj.seek(0)
        sep_to_use = manual_sep if (manual_sep and manual_sep != "Auto") else detect_csv_separator(file_obj)
        
        try:
            df = pd.read_csv(
                file_obj,
                header=header_row,
                nrows=nrows,
                sep=sep_to_use,
                engine='python',
                encoding='utf-8',
                on_bad_lines='skip'
            )
        except UnicodeDecodeError:
            file_obj.seek(0)
            df = pd.read_csv(
                file_obj,
                header=header_row,
                nrows=nrows,
                sep=sep_to_use,
                engine='python',
                encoding='cp1250',
                on_bad_lines='skip'
            )
            
        if header_row is not None and skip_after > 0:
            df = df.iloc[int(skip_after):].reset_index(drop=True)
            
        return df
    return pd.DataFrame()

def parse_numeric_series(series: pd.Series) -> np.ndarray:
    """Konwertuje serię na float64, usuwając szum tekstowy i zamieniając przecinki dziesiętne."""
    if series is None or len(series) == 0:
        return np.array([], dtype=np.float64)
        
    if pd.api.types.is_numeric_dtype(series):
        return series.to_numpy(dtype=np.float64)
        
    s_str = series.astype(str).str.strip()
    s_str = s_str.str.replace(',', '.', regex=False)
    s_str = s_str.apply(lambda x: re.sub(r'[^\d\.\+\-eE]', '', x) if isinstance(x, str) else x)
    
    return pd.to_numeric(s_str, errors='coerce').to_numpy(dtype=np.float64)

# --- SILNIK OBLICZENIOWY ISO 6892-1 / ASTM E8M ---
def process_tensile_dataset(
    df_raw: pd.DataFrame,
    force_col: str,
    disp_col: str,
    force_unit: str,
    disp_unit: str,
    l0: float,
    s0: float,
    sample_id: str
):
    """Deterministyczny silnik wyznaczania parametrów wytrzymałościowych z kompensacją toe."""
    if force_col not in df_raw.columns or disp_col not in df_raw.columns:
        return None
        
    raw_force = parse_numeric_series(df_raw[force_col])
    raw_disp = parse_numeric_series(df_raw[disp_col])
    
    clean_mask = np.isfinite(raw_force) & np.isfinite(raw_disp)
    force_arr = raw_force[clean_mask]
    disp_arr = raw_disp[clean_mask]
    
    if len(force_arr) < 20:
        return None

    # Normalizacja do jednostek SI: Siła [N], Droga [mm]
    if force_unit == "kN":
        force_arr = force_arr * 1e3
    elif force_unit == "MN":
        force_arr = force_arr * 1e6
        
    if disp_unit in ["μm", "um"]:
        disp_arr = disp_arr / 1e3
    elif disp_unit == "m":
        disp_arr = disp_arr * 1e3

    stress = force_arr / s0
    strain = disp_arr / l0
    
    rm_idx = int(np.argmax(stress))
    if rm_idx < 5:
        rm_idx = len(stress) - 1
        
    fm_kn = float(force_arr[rm_idx] / 1000.0)
    rm_mpa = float(stress[rm_idx])
    fu_kn = float(force_arr[-1] / 1000.0)
    ru_mpa = float(stress[-1])
    a_pct = float(strain[-1] * 100.0)
    
    # Pasmo odcięcia szumów brzegowych i uślizgu uchwytów (15% - 75% Rm)
    stress_min_cutoff = 0.15 * rm_mpa
    stress_max_cutoff = 0.75 * rm_mpa
    
    search_indices = np.where(
        (stress >= stress_min_cutoff) & 
        (stress <= stress_max_cutoff) & 
        (np.arange(len(stress)) <= rm_idx)
    )[0]
    
    if len(search_indices) < 15:
        search_indices = np.arange(0, max(15, int(rm_idx * 0.75)))

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

    # Fallback dla silnie zaszumionych sygnałów
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

    if best_slope <= 0.0:
        delta_e = strain[rm_idx] - strain[0]
        best_slope = (stress[rm_idx] - stress[0]) / delta_e if delta_e > 0 else 1000.0
        best_intercept = stress[0] - best_slope * strain[0]

    # Kompensacja luzu początkowego (Toe compensation)
    toe_offset = -best_intercept / best_slope
    strain_corr = strain - toe_offset
    
    # Wyznaczenie granicy plastyczności Rp0.2 (offset 0.2%)
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
        "sample_id": str(sample_id),
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
        "toe_offset_pct": toe_offset * 100.0,
        "valid_points_count": len(force_arr)
    }

# --- INTERFEJS UŻYTKOWNIKA STREAMLIT ---
st.title("🔬 Moduł Analizy Próby Rozciągania")
st.caption("Standardy metrologiczne: PN-EN ISO 6892-1 / ASTM E8M (Obsługa Excel & CSV z Kompensacją Luzu)")

with st.sidebar:
    st.header("⚙️ 1. Pliki Wejściowe")
    uploaded_files = st.file_uploader(
        "Załaduj pliki pomiarowe (.xlsx, .xls, .csv)",
        type=["xlsx", "xls", "csv"],
        accept_multiple_files=True
    )

if not uploaded_files:
    st.info("👆 Przeciągnij plik(i) Excel lub CSV w lewym panelu, aby rozpocząć analizę.")
else:
    has_csv = any(f.name.lower().endswith('.csv') for f in uploaded_files)
    
    sample_registry = {}
    for f in uploaded_files:
        if f.name.lower().endswith(('.xlsx', '.xls')):
            xl = pd.ExcelFile(f)
            for sheet in xl.sheet_names:
                sample_registry[f"{f.name} -> {sheet}"] = {"file": f, "item": sheet, "type": "excel"}
        else:
            sample_registry[f.name] = {"file": f, "item": f.name, "type": "csv"}
            
    sample_keys = list(sample_registry.keys())

    with st.sidebar:
        st.markdown("---")
        st.header("🛠️ 2. Struktura Danych")
        
        active_sep = "Auto"
        if has_csv:
            csv_sep_choice = st.selectbox(
                "Separator kolumn dla CSV:",
                ["Auto", "; (Średnik - PL/DE)", ", (Przecinek - EN/US)", "\\t (Tabulator)"],
                index=0,
                help="Wybierz separator używany w Twoim pliku CSV."
            )
            sep_map = {
                "Auto": "Auto",
                "; (Średnik - PL/DE)": ";",
                ", (Przecinek - EN/US)": ",",
                "\\t (Tabulator)": "\t"
            }
            active_sep = sep_map[csv_sep_choice]
        
        header_row_idx = st.number_input(
            "Wiersz z nagłówkami (0 = 1. wiersz):",
            value=5,
            min_value=0,
            max_value=50,
            step=1
        )
        skip_after_header = st.number_input(
            "Pomiń wiersze pod nagłówkiem (np. jednostki):",
            value=1,
            min_value=0,
            max_value=10,
            step=1
        )

    # SEKCJA DIAGNOSTYCZNA / RAW INSPECTOR
    with st.expander("🔍 Podgląd Surowego Nagłówka i Metadanych (Wskaźnik Wierszy)", expanded=True):
        col_insp1, col_insp2 = st.columns([3, 1])
        with col_insp1:
            selected_preview_key = st.selectbox("Wybierz próbkę do inspekcji:", sample_keys)
        with col_insp2:
            n_rows_preview = st.number_input("Wiersze do podglądu:", min_value=3, max_value=30, value=10, step=1)
            
        target_meta = sample_registry[selected_preview_key]
        df_raw_preview = load_raw_dataframe(
            target_meta["file"], target_meta["item"], header_row=None, nrows=int(n_rows_preview), manual_sep=active_sep
        )
        df_raw_preview.index = [f"Wiersz {i}" for i in range(len(df_raw_preview))]
        st.dataframe(df_raw_preview, use_container_width=True)

    with st.sidebar:
        st.markdown("---")
        st.header("📐 3. Geometria Próbki")
        col_g1, col_g2 = st.columns(2)
        with col_g1:
            manual_l0 = st.number_input("Baza L₀ [mm]", value=40.0, step=1.0, min_value=1.0)
        with col_g2:
            manual_s0 = st.number_input("Przekrój S₀ [mm²]", value=73.125, step=0.5, min_value=0.1)

    # Odczyt nagłówków próbki referencyjnej
    df_sample_preview = load_raw_dataframe(
        target_meta["file"], target_meta["item"], 
        header_row=int(header_row_idx), skip_after=int(skip_after_header), manual_sep=active_sep
    )
    df_sample_preview.columns = [str(col).strip() for col in df_sample_preview.columns]
    column_options = df_sample_preview.columns.tolist()

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
            "Wybierz zakres raportu:",
            ["Pojedyncza Próbka (Raport Szczegółowy)", "Zbiorcza Seria Pomiarowa (Multi-Sample)"]
        )

    # =========================================================================
    # TRYB 1: POJEDYNCZA PRÓBKA (RAPORT SZCZEGÓŁOWY)
    # =========================================================================
    if analysis_mode == "Pojedyncza Próbka (Raport Szczegółowy)":
        with st.sidebar:
            selected_single_key = st.selectbox("Wybierz próbkę do analizy:", sample_keys)
            
        t_meta = sample_registry[selected_single_key]
        df_target = load_raw_dataframe(
            t_meta["file"], t_meta["item"], 
            header_row=int(header_row_idx), skip_after=int(skip_after_header), manual_sep=active_sep
        )
        df_target.columns = [str(col).strip() for col in df_target.columns]
        
        res = process_tensile_dataset(
            df_target, force_column, disp_column, force_unit, disp_unit, manual_l0, manual_s0, selected_single_key
        )
        
        if res is None:
            st.error(
                f"Błąd: Nie udało się wyodrębnić punktów pomiarowych z kolumn '{force_column}' i '{disp_column}'. "
                "Upewnij się w lewym panelu, że wskaźnik wiersza nagłówka oraz separator są ustawione poprawnie."
            )
        else:
            col1, col2, col3, col4, col5 = st.columns(5)
            col1.metric("Wytrzymałość Doraźna Rm", f"{res['Rm_MPa']:.1f} MPa")
            col2.metric("Granica Plastyczności Rp0.2", f"{res['Rp02_MPa']:.1f} MPa" if res['Rp02_MPa'] else "Brak")
            col3.metric("Naprężenie Zerwania Ru", f"{res['Ru_MPa']:.1f} MPa")
            col4.metric("Wydłużenie Całkowite A", f"{res['A_pct']:.1f} %")
            col5.metric("Maksymalna Siła Fm", f"{res['Fm_kN']:.2f} kN")
            
            # Anchor krzywej do początku układu (0,0) – odcięcie wartości ujemnych po kompensacji
            plot_mask = res["strain_corr_pct"] >= 0
            plot_strain = np.insert(res["strain_corr_pct"][plot_mask], 0, 0.0)
            plot_stress = np.insert(res["stress"][plot_mask], 0, 0.0)

            fig = go.Figure()
            
            # 1. Główna krzywa rozciągania
            fig.add_trace(go.Scatter(
                x=plot_strain,
                y=plot_stress,
                mode='lines',
                name=f'Krzywa σ-ε ({selected_single_key})',
                line=dict(color='#0066cc', width=2.5),
                cliponaxis=False
            ))
            
            # 2. Ograniczona prosta offsetowa 0.2%
            rp02_target = res["Rp02_MPa"] if res["Rp02_MPa"] else (res["Rm_MPa"] * 0.8)
            max_offset_stress = min(res["Rm_MPa"] * 0.95, rp02_target * 1.15)
            max_offset_strain = (max_offset_stress / res["slope"])
            
            e_range = np.linspace(0, max_offset_strain, 50)
            fig.add_trace(go.Scatter(
                x=(e_range + 0.002) * 100,
                y=res["slope"] * e_range,
                mode='lines',
                name='Offset 0.2% (Sztywność układu)',
                line=dict(color='rgba(230, 0, 0, 0.75)', dash='dash', width=1.8),
                cliponaxis=False
            ))
            
            # 3. Punkt Rp0.2 (bez dublowania w legendzie)
            if res["Rp02_MPa"] and res["Rp02_strain_pct"]:
                fig.add_trace(go.Scatter(
                    x=[res["Rp02_strain_pct"]],
                    y=[res["Rp02_MPa"]],
                    mode='markers+text',
                    text=[f"<b>Rp0.2:</b> {res['Rp02_MPa']:.1f} MPa"],
                    textposition="top left",
                    marker=dict(color='#d90429', size=9, symbol='diamond'),
                    showlegend=False,
                    cliponaxis=False
                ))
                
            # 4. Punkt Rm (bez dublowania w legendzie)
            fig.add_trace(go.Scatter(
                x=[res["Rm_strain_pct"]],
                y=[res["Rm_MPa"]],
                mode='markers+text',
                text=[f"<b>Rm:</b> {res['Rm_MPa']:.1f} MPa"],
                textposition="top center",
                marker=dict(color='#f77f00', size=10, symbol='circle'),
                showlegend=False,
                cliponaxis=False
            ))
            
            # 5. Punkt Ru (bez dublowania w legendzie)
            fig.add_trace(go.Scatter(
                x=[res["strain_corr_pct"][-1]],
                y=[res["Ru_MPa"]],
                mode='markers+text',
                text=[f"<b>Ru:</b> {res['Ru_MPa']:.1f} MPa (A = {res['A_pct']:.1f}%)"],
                textposition="bottom left",
                marker=dict(color='#8d99ae', size=9, symbol='x'),
                showlegend=False,
                cliponaxis=False
            ))
            
            # Limity osi: start bezwzględny od 0, headroom +15% na Y zapobiegający ucinaniu etykiet
            y_headroom = float(res["Rm_MPa"] * 1.15)
            x_headroom = float(res["strain_corr_pct"][-1] * 1.06)
            
            fig.update_layout(
                title=dict(
                    text=f"<b>Charakterystyka Rozciągania – {selected_single_key}</b>", 
                    font=dict(size=17, color="#f8f9fa")
                ),
                xaxis=dict(
                    title="Odkształcenie skorygowane ε_corr [%]",
                    range=[0, x_headroom],
                    gridcolor="rgba(255, 255, 255, 0.1)",
                    zeroline=True,
                    zerolinecolor="rgba(255, 255, 255, 0.2)"
                ),
                yaxis=dict(
                    title="Naprężenie inżynierskie σ [MPa]",
                    range=[0, y_headroom],
                    gridcolor="rgba(255, 255, 255, 0.1)",
                    zeroline=True,
                    zerolinecolor="rgba(255, 255, 255, 0.2)"
                ),
                hovermode="closest",
                template="plotly_dark",
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=1.02,
                    xanchor="right",
                    x=1.0,
                    bgcolor="rgba(0,0,0,0)"
                ),
                height=580,
                margin=dict(l=60, r=40, t=70, b=50)
            )
            
            st.plotly_chart(fig, use_container_width=True)

    # =========================================================================
    # TRYB 2: ZBIORCZA SERIA POMIAROWA (MULTI-SAMPLE)
    # =========================================================================
    else:
        with st.sidebar:
            selected_multi_keys = st.multiselect(
                "Wybierz próbki do zestawienia serii:",
                sample_keys,
                default=sample_keys
            )
            
        if not selected_multi_keys:
            st.warning("Wybierz co najmniej jedną próbkę do wygenerowania raportu serii.")
        else:
            batch_results = []
            fig_cum = go.Figure()
            max_rm_overall = 0.0
            
            for key in selected_multi_keys:
                t_meta = sample_registry[key]
                df_item = load_raw_dataframe(
                    t_meta["file"], t_meta["item"], 
                    header_row=int(header_row_idx), skip_after=int(skip_after_header), manual_sep=active_sep
                )
                df_item.columns = [str(col).strip() for col in df_item.columns]
                
                res_item = process_tensile_dataset(
                    df_item, force_column, disp_column, force_unit, disp_unit, manual_l0, manual_s0, key
                )
                
                if res_item is not None:
                    if res_item["Rm_MPa"] > max_rm_overall:
                        max_rm_overall = res_item["Rm_MPa"]
                        
                    batch_results.append({
                        "Próbka": key,
                        "Fm [kN]": round(res_item["Fm_kN"], 2),
                        "Rm [MPa]": round(res_item["Rm_MPa"], 1),
                        "Rp0.2 [MPa]": round(res_item["Rp02_MPa"], 1) if res_item["Rp02_MPa"] else np.nan,
                        "Fu [kN]": round(res_item["Fu_kN"], 2),
                        "Ru [MPa]": round(res_item["Ru_MPa"], 1),
                        "A [%]": round(res_item["A_pct"], 1)
                    })
                    
                    b_mask = res_item["strain_corr_pct"] >= 0
                    b_strain = np.insert(res_item["strain_corr_pct"][b_mask], 0, 0.0)
                    b_stress = np.insert(res_item["stress"][b_mask], 0, 0.0)
                    
                    fig_cum.add_trace(go.Scatter(
                        x=b_strain,
                        y=b_stress,
                        mode='lines',
                        name=key,
                        line=dict(width=1.8),
                        cliponaxis=False
                    ))
                    
            if not batch_results:
                st.error("Żaden ze wskazanych plików nie zawierał poprawnych danych do analizy.")
            else:
                df_batch = pd.DataFrame(batch_results)
                
                st.subheader("📈 Skumulowany Przebieg Krzywych Rozciągania Serii")
                fig_cum.update_layout(
                    xaxis=dict(
                        title="Odkształcenie skorygowane ε_corr [%]",
                        range=[0, None],
                        gridcolor="rgba(255, 255, 255, 0.1)"
                    ),
                    yaxis=dict(
                        title="Naprężenie inżynierskie σ [MPa]",
                        range=[0, float(max_rm_overall * 1.15)],
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