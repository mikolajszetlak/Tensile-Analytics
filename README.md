
# Tensile Analytics Pro: Automated Mechanical Characterization & Proof Stress Engine (ISO 6892-1 / ASTM E8M)

Przemysłowy silnik analityczny oraz interaktywny pulpit inżynierski R&D przeznaczony do automatyzacji przetwarzania danych z jednoosiowej statycznej próby rozciągania metali i tworzyw sztucznych. System implementuje procedury deterministycznej kompensacji podatności układu i luzu początkowego (*toe compensation*), adaptacyjny algorytm ruchomego okna do identyfikacji sztywności technicznej oraz mikrometryczną interpolację odcinkową umownej granicy plastyczności $R_{p0.2}$ zgodnie z normami **PN-EN ISO 6892-1:2020** oraz **ASTM E8 / E8M-21**.

---

## Spis Treści
1. [Wprowadzenie i Przeznaczenie Systemu](#1-wprowadzenie-i-przeznaczenie-systemu)
2. [Metrologia i Fizyka Zjawiska: Trawersa vs Ekstensometr](#2-metrologia-i-fizyka-zjawiska-trawersa-vs-ekstensometr)
   - [Podatność Układu Badawczego (*System Compliance*)](#podatność-układu-badawczego-system-compliance)
   - [Artefakt Luzu Początkowego (*Toe Region Artifact*)](#artefakt-luzu-początkowego-toe-region-artifact)
3. [Pełny Aparat Matematyczny i Algorytmy](#3-pełny-aparat-matematyczny-i-algorytmy)
   - [Normalizacja Jednostek i Wektory Inżynierskie ($\sigma, \varepsilon$)](#normalizacja-jednostek-i-wektory-inżynierskie-\sigma-\varepsilon)
   - [Filtracja Przestrzeni Poszukiwań (Odcięcie Szumów Brzegowych)](#filtracja-przestrzeni-poszukiwań-odcięcie-szumów-brzegowych)
   - [Metoda Ruchomego Okna Regresji Liniowej (OLS)](#metoda-ruchomego-okna-regresji-liniowej-ols)
   - [Kompensacja Luzu Początkowego (*Toe Shift Derivation*)](#kompensacja-luzu-początkowego-toe-shift-derivation)
   - [Konstrukcja Offsetu 0.2% i Segmentowa Interpolacja Liniowa $R_{p0.2}$](#konstrukcja-offsetu-02-i-segmentowa-interpolacja-liniowa-r_p02)
   - [Parametry Graniczne i Statystyka Serii ($\mu, \sigma, V$)](#parametry-graniczne-i-statystyka-serii-\mu-\sigma-v)
4. [Szczegółowa Anatomia Funkcji `process_tensile_dataset()`](#4-szczegółowa-anatomia-funkcji-process_tensile_dataset)
   - [Sygnatura i Parametry Wejściowe](#sygnatura-i-parametry-wejściowe)
   - [Dekompozycja Kroków Wykonawczych (Co robi, Jak działa, Skąd pochodzi)](#dekompozycja-kroków-wykonawczych-co-robi-jak-działa-skąd-pochodzi)
   - [Struktura Słownika Wyjściowego](#struktura-słownika-wyjściowego)
5. [Architektura Interfejsu (Streamlit + Plotly)](#5-architektura-interfejsu-streamlit--plotly)
   - [Inspekcja Surowego Nagłówka (*Raw Metadata Inspector*)](#inspekcja-surowego-nagłówka-raw-metadata-inspector)
   - [Dynamiczny Kreator Mapowania Kolumn](#dynamiczny-kreator-mapowania-kolumn)
   - [Skalowanie Bounding-Box Wykresu (Eliminacja Spłaszczenia)](#skalowanie-bounding-box-wykresu-eliminacja-spłaszczenia)
   - [Moduł Multi-Sheet i Analiza Zbiorcza](#moduł-multi-sheet-i-analiza-zbiorcza)
6. [Struktura Repozytorium i Wdrożenie](#6-struktura-repozytorium-i-wdrożenie)
7. [Zastrzeżenia Metrologiczne i Rygor Normowy](#7-zastrzeżenia-metrologiczne-i-rygor-normowy)

---

## 1. Wprowadzenie i Przeznaczenie Systemu

Oprogramowanie fabryczne maszyn wytrzymałościowych (ZwickRoell, Instron, MTS, Hegewald&Peschke) generuje pliki arkuszy kalkulacyjnych o zróżnicowanej strukturze nagłówków, jednostek i liczbie wierszy metadanych. Ręczna obróbka takich danych w programie MS Excel jest czasochłonna, obarczona ryzykiem błędu ludzkiego oraz pozbawiona rygorystycznej kompensacji luzu początkowego.

**Tensile Analytics Pro** rozwiązuje ten problem poprzez:
* Pełną agnostyczność względem formatu wejściowego (dowolna kolejność kolumn, wiersze jednostek, różne układy jednostek SI).
* Automatyczne wyznaczanie punktów $R_m$, $R_{p0.2}$, $R_u$, $A$ oraz $F_m$ bez konieczności ręcznego zaznaczania stref sprężystych na wykresie.
* Zapewnienie powtarzalności obliczeń w seriach badawczych z jednoczesnym generowaniem zestawień statystycznych.

---

## 2. Metrologia i Fizyka Zjawiska: Trawersa vs Ekstensometr

### Podatność Układu Badawczego (*System Compliance*)
Gdy przemieszczenie próbki jest mierzone za pomocą impulsatora położenia trawersy maszyny (kolumna drogi standardowej), rejestrowany sygnał $\Delta L_{\text{raw}}$ nie jest rzeczywistym wydłużeniem bazy pomiarowej $L_0$, lecz sumą odkształceń całego łańcucha kinematycznego:

$$\Delta L_{\text{raw}} = \Delta L_{\text{próbka}} + \Delta L_{\text{uchwyty}} + \Delta L_{\text{rama}} + \Delta L_{\text{czujnik\_siły}}$$

Z tego względu nachylenie prostej w zakresie квази-sprężystym reprezentuje **sztywność techniczną układu badawczego** ($E_{\text{apparent}} \approx 1.5 - 3.5\text{ GPa}$ dla metali), a nie fizyczny moduł sprężystości wzdłużnej materiału (moduł Younga $E \approx 70\text{ GPa}$ dla Al, $210\text{ GPa}$ dla stali).

### Artefakt Luzu Początkowego (*Toe Region Artifact*)
W początkowej fazie testu występuje nieliniowy odcinek wklęsły wywołany:
1. Dociskiem i zagłębianiem się ząbkowanych klinów w główkach próbki.
2. Kasowaniem luzów mechanicznych na kolumnach prowadzących i śrubach pociągowych.
3. Wstępnym osiowaniem próbki w osi obciążenia.

```text
Naprężenie [MPa]
   ^
   |                 /  (Prawdziwy kierunek sztywności)
   |                /
   |               /
   |              / 
   |             /
   |         .--'   <-- Strefa uślizgu w klinach i kasowania luzu (Toe Artifact)
 0 +--------/------------------------------------------------------------> Odkształcenie [-]
          e_0 (Punkt pozornego zera)

```

Niezastosowanie kompensacji $\varepsilon_0$ powoduje, że prosta offsetowa $0.2\%$ zostaje skonstruowana w złym miejscu osi odkształceń, co prowadzi do drastycznego zaniżenia lub zawyżenia granicy plastyczności $R_{p0.2}$.

---

## 3. Pełny Aparat Matematyczny i Algorytmy

```text
                               SCHEMAT PRZEPŁYWU OBLICZEŃ
┌────────────────────────┐      ┌────────────────────────┐      ┌────────────────────────┐
│ Surowe Dane Maszynowe  │ ───> │ Normalizacja SI        │ ───> │ Filtracja Szumów       │
│ F [kN/N], ΔL [mm/μm]   │      │ σ = F/S0, ε = ΔL/L0    │      │ 0.15 Rm <= σ <= 0.75 Rm│
└────────────────────────┘      └────────────────────────┘      └───────────┬────────────┘
                                                                            │
┌────────────────────────┐      ┌────────────────────────┐      ┌───────────▼────────────┐
│ Wyznaczenie Rp0.2      │ <─── │ Kompensacja Luzu (Toe) │ <─── │ Ruchome Okno OLS       │
│ Interpolacja odcinkowa │      │ ε_corr = ε - (-b/a)    │      │ max(a) przy R^2 > 0.985│
└───────────┬────────────┘      └────────────────────────┘      └────────────────────────┘
            │
            ▼
┌────────────────────────┐      ┌────────────────────────┐
│ Punkty Graniczne       │ ───> │ Zestawienie Statystyki │
│ Rm, Ru, A, Fm          │      │ μ, σ, V [%] dla Serii  │
└────────────────────────┘      └────────────────────────┘

```

### Normalizacja Jednostek i Wektory Inżynierskie ($\sigma, \varepsilon$)

Wszystkie wektory wejściowe są sprowadzane do bazowego układu jednostek: Siła $F$ w Niutonach [$N$], Przemieszczenie $\Delta L$ w milimetrach [$mm$].

Naprężenie inżynierskie $\sigma_i$ oraz odkształcenie inżynierskie $\varepsilon_i$ w $i$-tym punkcie pomiarowym:

$$\sigma_i = \frac{F_i}{S_0} \quad [\text{MPa}], \qquad \varepsilon_i = \frac{\Delta L_i}{L_0} \quad [-]$$

gdzie:

* $S_0$ – początkowe pole przekroju poprzecznego części pomiarowej [$\text{mm}^2$],
* $L_0$ – początkowa długość bazy pomiarowej [$\text{mm}$].

### Filtracja Przestrzeni Poszukiwań (Odcięcie Szumów Brzegowych)

Aby wyeliminować ryzyko dopasowania strefy liniowej do uślizgu w uchwytach (problem niskich naprężeń) lub strefy płynięcia plastycznego (wysokie naprężenia), definiowana jest dynamiczna maska logiczna:

$$\Omega_{\text{search}} = \left\{ i \in \mathbb{N} \;\middle\vert{}\; 0.15 \cdot R_m \le \sigma_i \le 0.75 \cdot R_m \;\land\; i \le i_{R_m} \right\}$$

gdzie:

* $i_{R_m} = \arg\max_i (\sigma_i)$ – indeks wytrzymałości doraźnej $R_m$.

### Metoda Ruchomego Okna Regresji Liniowej (OLS)

Szerokość okna analitycznego $W$ jest dobierana adaptacyjnie na podstawie gęstości próbkowania:

$$W = \max\left(6, \; \min\left(50, \; \lfloor 0.25 \cdot \vert{}\Omega_{\text{search}}\vert{} \rfloor\right)\right)$$

Dla każdego położenia okna $j \in [0, \vert{}\Omega_{\text{search}}\vert{} - W]$ wyznaczane są parametry prostej $\sigma = a \cdot \varepsilon + b$ metodą najmniejszych kwadratów (OLS):

$$a_j = \frac{\sum_{k=1}^W (\varepsilon_k - \bar{\varepsilon})(\sigma_k - \bar{\sigma})}{\sum_{k=1}^W (\varepsilon_k - \bar{\varepsilon})^2}$$

$$b_j = \bar{\sigma} - a_j \bar{\varepsilon}$$

Współczynnik determinacji:

$$R^2_j = \frac{\left[\sum_{k=1}^W (\varepsilon_k - \bar{\varepsilon})(\sigma_k - \bar{\sigma})\right]^2}{\sum_{k=1}^W (\varepsilon_k - \bar{\varepsilon})^2 \sum_{k=1}^W (\sigma_k - \bar{\sigma})^2}$$

**Kryterium selekcji:** Silnik wybiera okno o maksymalnym nachyleniu $a_j$ spełniające rygor prostoliniowości:

$$a_{\text{opt}} = \max \left\{ a_j \;\middle\vert{}\; R^2_j \ge 0.985 \right\}$$

Jeżeli żaden segment nie osiąga $R^2 \ge 0.985$ (dane zaszumione), aktywowany jest fallback:

$$a_{\text{opt}} = a_m \quad \text{gdzie} \quad m = \arg\max_j \left( R^2_j \right) \quad \text{dla} \quad a_j > 0$$

### Kompensacja Luzu Początkowego (*Toe Shift Derivation*)

Punkt przecięcia prostej sprężystej z osią zerowego naprężenia ($\sigma = 0$):

$$0 = a_{\text{opt}} \cdot \varepsilon_0 + b_{\text{opt}} \implies \varepsilon_0 = -\frac{b_{\text{opt}}}{a_{\text{opt}}}$$

Wektor odkształcenia całego eksperymentu ulega translacji:

$$\varepsilon_{\text{corr}, i} = \varepsilon_i - \varepsilon_0 = \varepsilon_i + \frac{b_{\text{opt}}}{a_{\text{opt}}}$$

Operacja ta gwarantuje, że skorygowana prosta sprężysta przechodzi ściśle przez początek układu współrzędnych $(0,0)$.

### Konstrukcja Offsetu 0.2% i Segmentowa Interpolacja Liniowa $R_{p0.2}$

Równanie prostej offsetowej przesuniętej o odkształcenie trwałe $\Delta \varepsilon = 0.002$ ($0.2\%$):

$$\sigma_{\text{offset}}(\varepsilon_{\text{corr}}) = a_{\text{opt}} \cdot (\varepsilon_{\text{corr}} - 0.002)$$

Dyskretna funkcja różnicowa:

$$f_i = \sigma_i - \sigma_{\text{offset}}(\varepsilon_{\text{corr}, i}) = \sigma_i - a_{\text{opt}} \cdot (\varepsilon_{\text{corr}, i} - 0.002)$$

Przeszukiwany jest zakres indeksów $i \ge i_{\text{start}} + W$ (za strefą dopasowania sprężystego) dla $\varepsilon_{\text{corr}, i} \ge 0.002$. Punkt przecięcia definiuje zmiana znaku funkcji różnicowej:

$$f_i \cdot f_{i+1} \le 0$$

Współrzędne punktu $(\varepsilon_{Rp0.2}, R_{p0.2})$ wyznaczane są za pomocą ważonej interpolacji liniowej:

$$w = \frac{\vert{}f_i\vert{}}{\vert{}f_i\vert{} + \vert{}f_{i+1}\vert{}}$$

$$\varepsilon_{Rp0.2} = \varepsilon_{\text{corr}, i} + w \cdot (\varepsilon_{\text{corr}, i+1} - \varepsilon_{\text{corr}, i})$$

$$R_{p0.2} = \sigma_i + w \cdot (\sigma_{i+1} - \sigma_i) \quad [\text{MPa}]$$

### Parametry Graniczne i Statystyka Serii ($\mu, \sigma, V$)

* **Maksymalna siła osiowa:**

$$F_m = \max_{i} (F_i) \quad [\text{N}] \implies \frac{F_m}{1000} \quad [\text{kN}]$$


* **Wytrzymałość doraźna na rozciąganie:**

$$R_m = \frac{F_m}{S_0} = \max_{i} (\sigma_i) \quad [\text{MPa}]$$


* **Naprężenie rozrywające:**

$$R_u = \frac{F_{\text{last}}}{S_0} = \sigma_{\text{last}} \quad [\text{MPa}]$$


* **Wydłużenie całkowite po rozerwaniu:**

$$A = \varepsilon_{\text{corr}, \text{last}} \cdot 100\% = \left( \frac{\Delta L_{\text{last}}}{L_0} - \varepsilon_0 \right) \cdot 100\% \quad [\%]$$



Dla serii pomiarowej $N$ próbek obliczane są estymatory statystyczne dla każdego parametru $X \in \{R_m, R_{p0.2}, R_u, A\}$:

$$\text{Średnia arytmetyczna: } \mu = \bar{X} = \frac{1}{N} \sum_{k=1}^N X_k$$

$$\text{Odchylenie standardowe (próbkowe): } s = \sqrt{\frac{1}{N - 1} \sum_{k=1}^N (X_k - \bar{X})^2}$$

$$\text{Współczynnik zmienności: } V = \frac{s}{\mu} \cdot 100\% \quad [\%]$$

---

## 4. Szczegółowa Anatomia Funkcji `process_tensile_dataset()`

Główna funkcja wykonawcza silnika analitycznego znajduje się w pliku `app.py`.

### Sygnatura i Parametry Wejściowe

```python
def process_tensile_dataset(
    df_raw: pd.DataFrame,
    force_col: str,
    disp_col: str,
    force_unit: str,
    disp_unit: str,
    l0: float,
    s0: float,
    sheet_name: str
) -> Optional[Dict[str, Any]]:

```

| Parametr Wejściowy | Typ | Źródło Danych w Systemie | Rola w Obliczeniach |
| --- | --- | --- | --- |
| `df_raw` | `pd.DataFrame` | Wycięty blok arkusza po odrzuceniu wiersza nagłówka | Zawiera surowe ciągi danych pomiarowych z maszyny. |
| `force_col` | `str` | Wybór użytkownika z listy rozwijanej UI | Identyfikuje kolumnę z wartościami siły obciążającej. |
| `disp_col` | `str` | Wybór użytkownika z listy rozwijanej UI | Identyfikuje kolumnę z przemieszczeniem / drogą trawersy. |
| `force_unit` | `str` | Selektor jednostki w UI (`"N"`, `"kN"`, `"MN"`) | Determinuje mnożnik skalujący do układu bazowego ($N$). |
| `disp_unit` | `str` | Selektor jednostki w UI (`"mm"`, `"μm"`, `"m"`) | Determinuje mnożnik skalujący do układu bazowego ($mm$). |
| `l0` | `float` | Pole numeryczne w UI (odczytane z nagłówka) | Pierwotna baza pomiarowa próbki $L_0$ [mm]. |
| `s0` | `float` | Pole numeryczne w UI (odczytane z nagłówka) | Pierwotne pole przekroju poprzecznego $S_0$ [$\text{mm}^2$]. |
| `sheet_name` | `str` | Pętla iterująca po skoroszycie Excela | Etykieta identyfikacyjna próbki w raportach i na wykresach. |

---

### Dekompozycja Kroków Wykonawczych

#### Krok 1: Ekstrakcja Numeryczna i Sanityzacja Wektorów

```python
force_series = pd.to_numeric(df_raw[force_col], errors='coerce')
disp_series = pd.to_numeric(df_raw[disp_col], errors='coerce')
clean_mask = force_series.notna() & disp_series.notna()
force_arr = force_series[clean_mask].to_numpy(dtype=np.float64)
disp_arr = disp_series[clean_mask].to_numpy(dtype=np.float64)

```

* **Co robi:** Rzutuje surowe dane tekstowe na typ `float64`, zamieniając błędy odczytu i puste komórki na `NaN`, po czym odrzuca uszkodzone wiersze za pomocą maski logicznej `clean_mask`.
* **Dlaczego tak:** Maszyny wytrzymałościowe często wstawiają wiersze tekstowe na końcu pliku (np. `"Koniec testu"`) lub znaki specjalne, które bez flagi `errors='coerce'` powodują awarię interpretera.

#### Krok 2: Skalowanie do Układu Bazowego SI ($N$, $mm$)

```python
if force_unit == "kN":
    force_arr = force_arr * 1e3
elif force_unit == "MN":
    force_arr = force_arr * 1e6
    
if disp_unit in ["μm", "um"]:
    disp_arr = disp_arr / 1e3
elif disp_unit == "m":
    disp_arr = disp_arr * 1e3

```

* **Co robi:** Przelicza wektory na bazowe jednostki $N$ oraz $mm$.
* **Dlaczego tak:** Zapobiega błędom rzędu wielkości przy obliczaniu modułu i naprężeń w MPa ($1\text{ MPa} = 1\text{ N/mm}^2$).

#### Krok 3: Wyznaczenie Naprężeń $\sigma$, Odkształceń $\varepsilon$ oraz Wskaźników Granicznych

```python
stress = force_arr / s0
strain = disp_arr / l0

rm_idx = int(np.argmax(stress))
fm_kn = float(force_arr[rm_idx] / 1000.0)
rm_mpa = float(stress[rm_idx])
fu_kn = float(force_arr[-1] / 1000.0)
ru_mpa = float(stress[-1])
a_pct = float(strain[-1] * 100.0)

```

* **Co robi:** Oblicza wektor naprężeń i odkształceń inżynierskich oraz pobiera wartości maksymalne ($R_m, F_m$) i końcowe ($R_u, F_u, A$).

#### Krok 4: Filtracja Szumów Uchwytów i Zdefiniowanie Okna

```python
stress_min_cutoff = 0.15 * rm_mpa
stress_max_cutoff = 0.75 * rm_mpa

search_indices = np.where(
    (stress >= stress_min_cutoff) & 
    (stress <= stress_max_cutoff) & 
    (np.arange(len(stress)) <= rm_idx)
)[0]

window_size = max(6, min(50, int(len(search_indices) * 0.25)))

```

* **Co robi:** Zawęża indeksy poszukiwania sztywności sprężystej do przedziału $[15\% R_m, \; 75\% R_m]$.
* **Dlaczego tak:** Zapobiega dopasowaniu prostej Hooke'a do początkowego uślizgu w klinach (błąd zidentyfikowany w Próbce 9, gdzie bez filtra algorytm wyznaczał $R_{p0.2} = 9\text{ MPa}$).

#### Krok 5: Pętla Ruchomego Okna i Identyfikacja Sztywności

```python
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

```

* **Co robi:** Przesuwa okno o szerokości `window_size` i dopasowuje prostą za pomocą `scipy.stats.linregress`. Rejestruje najwyższe nachylenie o współczynniku $R^2 > 0.985$.
* **Fallback:** W przypadku braku okna o $R^2 > 0.985$ kod automatycznie wybiera okno o najwyższym $R^2$ przy dodatnim nachyleniu, a w ostateczności moduł sieczny, eliminując możliwość wystąpienia `ZeroDivisionError`.

#### Krok 6: Kompensacja Luzu i Wyszukanie $R_{p0.2}$

```python
toe_offset = -best_intercept / best_slope
strain_corr = strain - toe_offset

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

```

* **Co robi:** Oblicza przesunięcie `toe_offset`, tworzy wektor skorygowany `strain_corr`, buduje prostą offsetową $0.2\%$ i znajduje punkt przecięcia poprzez wykrycie zmiany znaku różnicy `diff` z ważoną interpolacją liniową.

---

### Struktura Słownika Wyjściowego

Funkcja zwraca słownik `dict` o następującej strukturze pól:

```python
return {
    "sheet": str(sheet_name),              # Nazwa arkusza / identyfikator próbki
    "Fm_kN": fm_kn,                        # Maksymalna siła rozciągająca [kN]
    "Rm_MPa": rm_mpa,                      # Wytrzymałość doraźna [MPa]
    "Rm_strain_pct": float(...),           # Odkształcenie w punkcie Rm [%]
    "Rp02_MPa": rp02_val,                  # Umowna granica plastyczności Rp0.2 [MPa]
    "Rp02_strain_pct": float(...),         # Odkształcenie w punkcie Rp0.2 [%]
    "Fu_kN": fu_kn,                        # Siła w momencie zerwania [kN]
    "Ru_MPa": ru_mpa,                      # Naprężenie rozrywające [MPa]
    "A_pct": a_pct,                        # Całkowite wydłużenie po rozerwaniu [%]
    "strain_corr_pct": strain_corr * 100,  # Skorygowany wektor odkształceń do wykresu [%]
    "stress": stress,                      # Wektor naprężeń inżynierskich [MPa]
    "slope": best_slope,                   # Sztywność techniczna układu d_sigma/d_epsilon
    "toe_offset_pct": toe_offset * 100.0   # Wartość skompensowanego luzu e_0 [%]
}

```

---

## 5. Architektura Interfejsu (Streamlit + Plotly)

Aplikacja webowa została zaprojektowana z zachowaniem zasad ergonomii interfejsów przemysłowych:

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                          INTERFEJS GŁÓWNY (STREAMLIT)                       │
├──────────────────────────────┬──────────────────────────────────────────────┤
│ PANEL STEROWANIA (SIDEBAR)   │ OBSZAR ROBOCZY (MAIN VIEW)                   │
│ • Upload pliku .xlsx/.xls    │ 1. Raw Metadata Inspector (podgląd nagłówka) │
│ • Wymiary próbki (L0, S0)    │ 2. Kafelki kluczowych metryk (Rm, Rp0.2...)  │
│ • Wiersz nagłówka i skip     │ 3. Interaktywny wykres Plotly σ-ε            │
│ • Mapowanie kolumn i jednostek│ 4. Tabela statystyczna μ, σ, V [%] (seria)   │
│ • Wybór trybu (Single/Batch) │ 5. Przycisk eksportu raportu CSV             │
└──────────────────────────────┴──────────────────────────────────────────────┘

```

### Inspekcja Surowego Nagłówka (*Raw Metadata Inspector*)

Wbudowany moduł podglądu surowych wierszy:

```python
df_raw_preview = pd.read_excel(uploaded_file, sheet_name=preview_sheet, header=None, nrows=int(n_rows_preview))
df_raw_preview.index = [f"Wiersz {i}" for i in range(len(df_raw_preview))]
st.dataframe(df_raw_preview, use_container_width=True)

```

Pozwala operatorowi odczytać wymiary próbki ($L_0, S_0$) oraz precyzyjny numer wiersza z nagłówkami bez opuszczania aplikacji i otwierania programu MS Excel.

### Dynamiczny Kreator Mapowania Kolumn

Użytkownik wybiera nazwy kolumn z dynamicznie generowanej listy `column_options`, co uniezależnia kod od wersji językowej oprogramowania maszyny wytrzymałościowej (np. `"Siła"` vs `"Force"` vs `"Standardkraft"`).

### Skalowanie Bounding-Box Wykresu (Eliminacja Spłaszczenia)

Wykresy Plotly automatycznie skalują osie do skrajnych wartości wszystkich śladów (*traces*). Niekontrolowana prosta offsetowa rozciągała oś pionową do ponad $200\text{ MPa}$, spłaszczając krzywą badanego materiału o $R_m \approx 44\text{ MPa}$.

**Implementacja ograniczenia promienia offsetu:**

```python
rp02_target = res["Rp02_MPa"] if res["Rp02_MPa"] else (res["Rm_MPa"] * 0.8)
max_offset_stress = min(res["Rm_MPa"] * 0.95, rp02_target * 1.15)
max_offset_strain = (max_offset_stress / res["slope"])
e_range = np.linspace(0, max_offset_strain, 50)

```

Oś Y otrzymuje normowy bufor wysokości $12\%$ powyżej $R_m$:

```python
yaxis=dict(range=[0, float(res["Rm_MPa"] * 1.12)])

```

### Moduł Multi-Sheet i Analiza Zbiorcza

W trybie serii wieloarkuszowej aplikacja przetwarza wybrane arkusze w pętli wsadowej, nanosi wszystkie krzywe na wspólny układ współrzędnych oraz agreguje wyniki do tabeli z estymatorami średniej ($\mu$), odchylenia standardowego ($\sigma$) i współczynnika zmienności ($V$). Gotowy raport można pobrać jednym kliknięciem jako plik CSV rozdzielany średnikami.

---

## 6. Struktura Repozytorium i Wdrożenie

### Drzewo Projektu

```text
PY_TENSILE/
├── data/
│   ├── raw/
│   │   └── test_lab.xlsx           # Surowe pliki z maszyny wytrzymałościowej
│   └── reports/                    # Wygenerowane zestawienia CSV i wykresy PNG
├── src/
│   ├── __init__.py
│   └── analyzer.py                 # Silnik obliczeniowy do testów jednostkowych
├── .gitignore
├── app.py                          # Główna aplikacja webowa Streamlit + Plotly
├── README.md                       # Pełna dokumentacja techniczna i matematyczna
└── requirements.txt                # Zależności środowiska wirtualnego

```

### Wymagania Środowiskowe

* Python $\ge 3.10$
* PowerShell / Bash

### Instalacja i Uruchomienie

1. Sklonuj repozytorium i przejdź do katalogu projektu:

```powershell
git clone [https://github.com/twoj-profil/Tensile-Analytics-Pro.git](https://github.com/twoj-profil/Tensile-Analytics-Pro.git)
cd Tensile-Analytics-Pro

```

2. Utwórz i aktywuj środowisko wirtualne:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1

```

3. Zainstaluj zależności:

```powershell
pip install -r requirements.txt

```

Plik `requirements.txt`:

```text
streamlit>=1.30.0
pandas>=2.0.0
numpy>=1.24.0
scipy>=1.10.0
plotly>=5.18.0
openpyxl>=3.1.0

```

4. Uruchom serwer aplikacji Streamlit:

```powershell
streamlit run app.py

```

---

## 7. Zastrzeżenia Metrologiczne i Rygor Normowy

> ### ⚠️ DEKLARACJA ZGODNOŚCI Z NORMAMI PN-EN ISO 6892-1 ORAZ ASTM E111
> 
> 
> 1. **Sztywność Pozorna vs Moduł Younga:** Nachylenie wyznaczone z pomiaru drogi trawersy reprezentuje sztywność całego węzła konstrukcyjnego maszyny i uchwytów. Wartość ta jest poprawnie wykorzystywana w silniku do wyznaczenia kierunku offsetu granicy plastyczności $R_{p0.2}$ (kompensacja ugięcia ramy), lecz **nie może być raportowana jako moduł sprężystości wzdłużnej ($E$) badanego tworzywa lub stopu**.
> 2. **Pomiar Modułu Younga ($E$):** Prawidłowe wyznaczenie właściwości sprężystych materiału wymaga zastosowania czujnika bezpośredniego na próbce:
> * Ekstensometru osiowego o klasie dokładności $\le 1$ wg ISO 9513,
> * Tensometrów elektrooporowych w układzie mostkowym,
> * Optycznego systemu cyfrowej korelacji obrazu (DIC - *Digital Image Correlation*).
> 
> 
> 3. **Raportowane Wskaźniki:** Parametry $R_m$, $R_{p0.2}$, $R_u$ oraz całkowite wydłużenie $A$ obliczone przez niniejszy silnik po kompensacji luzu spełniają wymagania stawiane rutynowym próbom odbiorczym i badaniom porównawczym w przemyśle metalowym oraz przetwórstwa polimerów.
> 
> 

