import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

# 1. Konfiguracja geometryczna (zgodna z arkuszem lab)
L0 = 40.0        # Baza pomiarowa [mm]
S0 = 73.125      # Pierwotne pole przekroju poprzecznego [mm^2]
file_path = "data/raw/lab_tensile_test.xlsx"

if not os.path.exists(file_path):
    raise FileNotFoundError(f"Nie odnaleziono pliku: {file_path}")

xl = pd.ExcelFile(file_path)
results_list = []

plt.figure(figsize=(10, 6))

for sheet in xl.sheet_names:
    # Wczytanie danych z pominięciem nagłówka tekstowego
    df = pd.read_excel(file_path, sheet_name=sheet, header=5)
    df = df.iloc[1:].reset_index(drop=True)
    
    force_n = pd.to_numeric(df.iloc[:, 1], errors='coerce').dropna().to_numpy()
    disp_mm = pd.to_numeric(df.iloc[:, 2], errors='coerce').dropna().to_numpy()
    
    stress = force_n / S0
    strain = disp_mm / L0
    
    # Parametry graniczne
    rm_idx = int(np.argmax(stress))
    fm_kn = float(force_n[rm_idx] / 1000.0)
    rm_mpa = float(stress[rm_idx])
    
    fu_kn = float(force_n[-1] / 1000.0)
    ru_mpa = float(stress[-1])
    a_pct = float(strain[-1] * 100.0)
    
    # Wyznaczenie sztywności technicznej układu do offsetu Rp0.2
    window_size = max(20, min(100, int(rm_idx * 0.12)))
    best_slope = 0.0
    best_intercept = 0.0
    
    for i in range(0, int(rm_idx * 0.85) - window_size, 2):
        w_strain = strain[i : i + window_size]
        w_stress = stress[i : i + window_size]
        slope, intercept, r_val, _, _ = linregress(w_strain, w_stress)
        if (r_val ** 2) > 0.990 and slope > best_slope:
            best_slope = slope
            best_intercept = intercept
            
    # Kompensacja luzu i wyznaczenie przecięcia Rp0.2
    toe_offset = -best_intercept / best_slope
    strain_corr = strain - toe_offset
    
    sigma_offset = best_slope * (strain_corr - 0.002)
    diff = stress - sigma_offset
    
    rp02_val = np.nan
    valid_range = np.where(strain_corr >= 0.002)[0]
    for idx in valid_range[:-1]:
        if diff[idx] * diff[idx + 1] <= 0:
            weight = abs(diff[idx]) / (abs(diff[idx]) + abs(diff[idx + 1]))
            rp02_val = float(stress[idx] + weight * (stress[idx + 1] - stress[idx]))
            break
            
    results_list.append({
        "Próbka": f"Próbka {sheet}",
        "F_m [kN]": round(fm_kn, 2),
        "R_m [MPa]": round(rm_mpa, 1),
        "R_p0.2 [MPa]": round(rp02_val, 1),
        "F_u [kN]": round(fu_kn, 2),
        "R_u [MPa]": round(ru_mpa, 1),
        "A [%]": round(a_pct, 1)
    })
    
    plt.plot(strain_corr * 100, stress, lw=1.2, alpha=0.8, label=f"Próbka {sheet}")

# Zestawienie tabelaryczne serii pomiarowej
df_report = pd.DataFrame(results_list)

print("=" * 70)
print("RAPORT WYNIKÓW PRÓBY ROZCIĄGANIA (ISO 6892-1 / DANE Z TRAWERSY)")
print("=" * 70)
print(df_report.to_string(index=False))
print("-" * 70)

# Analiza statystyczna serii
stats_summary = pd.DataFrame({
    "Wskaźnik": ["Średnia (μ)", "Odchylenie std (σ)", "Współczynnik zmienności (V) [%]"],
    "R_m [MPa]": [
        f"{df_report['R_m [MPa]'].mean():.2f}",
        f"{df_report['R_m [MPa]'].std():.2f}",
        f"{(df_report['R_m [MPa]'].std() / df_report['R_m [MPa]'].mean()) * 100:.2f}%"
    ],
    "R_p0.2 [MPa]": [
        f"{df_report['R_p0.2 [MPa]'].mean():.2f}",
        f"{df_report['R_p0.2 [MPa]'].std():.2f}",
        f"{(df_report['R_p0.2 [MPa]'].std() / df_report['R_p0.2 [MPa]'].mean()) * 100:.2f}%"
    ],
    "A [%]": [
        f"{df_report['A [%]'].mean():.2f}",
        f"{df_report['A [%]'].std():.2f}",
        f"{(df_report['A [%]'].std() / df_report['A [%]'].mean()) * 100:.2f}%"
    ]
})
print(stats_summary.to_string(index=False))
print("=" * 70)

# Formatowanie wykresu zbiorczego
plt.title(r"Krzywe rozciągania serii próbek $\sigma - \epsilon_{corr}$", fontsize=12, fontweight='bold')
plt.xlabel(r"Odkształcenie skorygowane $\epsilon_{corr}$ [%]", fontsize=10)
plt.ylabel(r"Naprężenie inżynierskie $\sigma$ [MPa]", fontsize=10)
plt.xlim(0, 30)
plt.ylim(0, df_report["R_m [MPa]"].max() * 1.15)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()

os.makedirs("data/reports", exist_ok=True)
plt.savefig("data/reports/tensile_series_summary.png", dpi=200)
df_report.to_csv("data/reports/tensile_report_summary.csv", index=False, sep=";")
print("Zapisano raport CSV oraz wykres w folderze: data/reports/")
plt.show()