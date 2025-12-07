# Sistema di Simulazione Radar Avanzato

Creerò un sistema completo di simulazione radar che considera tutti gli aspetti fisici e geografici richiesti. Prima spiegherò le equazioni fondamentali, poi implementerò il sistema in Go.

## 1. Equazioni Fondamentali

### 1.1 Equazione Radar di Base
L'**equazione radar** calcola la potenza ricevuta da un bersaglio:

```
Pr = (Pt × Gt × Gr × λ² × σ) / ((4π)³ × R⁴ × L)
```

Dove:
- **Pr**: potenza ricevuta (W)
- **Pt**: potenza trasmessa (W)
- **Gt**: guadagno antenna trasmittente
- **Gr**: guadagno antenna ricevente
- **λ**: lunghezza d'onda (m) = c/f
- **σ**: sezione radar equivalente (RCS) del bersaglio (m²)
- **R**: distanza radar-bersaglio (m)
- **L**: fattore di perdita totale

### 1.2 Orizzonte Geometrico e Orizzonte Radar

**Orizzonte geometrico** (linea di vista):
```
d = √(2 × R_terra × h)
```

**Orizzonte radar** (considera rifrazione atmosferica):
```
d_radar = √(2 × K × R_terra × h)
```

Dove:
- **R_terra**: raggio terrestre ≈ 6371 km
- **h**: altezza antenna (m)
- **K**: fattore di rifrazione atmosferica ≈ 4/3 (1.33)

Per due punti a diverse altezze:
```
d_totale = √(2 × K × R_terra × h1) + √(2 × K × R_terra × h2)
```

### 1.3 Attenuazione Atmosferica

L'attenuazione dipende da temperatura, umidità e frequenza:

**Attenuazione totale**:
```
L_atm = L_gas + L_pioggia
```

**Attenuazione per gas (approssimazione ITU-R)**:
```
γ_gas = (7.19×10⁻³ × f² × P) / ((f² + 0.34)×(f² + 0.27×P²)) + 
        (0.028×ρ×f²) / ((f-22.2)² + 8.5 + (f-183.3)² + 9)
```

Dove:
- **f**: frequenza (GHz)
- **P**: pressione atmosferica (hPa)
- **ρ**: densità vapor d'acqua (g/m³)

### 1.4 Conversioni Geodetiche

**Da coordinate geografiche a distanza**:
```
Δlat = lat2 - lat1
Δlon = lon2 - lon1

a = sin²(Δlat/2) + cos(lat1) × cos(lat2) × sin²(Δlon/2)
c = 2 × atan2(√a, √(1-a))
d = R_terra × c
```

**Azimuth** (direzione):
```
y = sin(Δlon) × cos(lat2)
x = cos(lat1) × sin(lat2) - sin(lat1) × cos(lat2) × cos(Δlon)
azimuth = atan2(y, x)
```

### 1.5 Angolo di Elevazione

```
α = atan2(h_target - h_radar, d_horizontal)
```

## 2. Algoritmi Implementati

### A. Calcolo Visibilità Line-of-Sight (LOS)
1. Calcola distanza orizzontale tra radar e target
2. Calcola orizzonte radar per entrambi i punti
3. Verifica se distanza < orizzonte totale
4. Considera curvatura terrestre e rifrazione

### B. Calcolo Attenuazione Atmosferica
1. Determina parametri atmosferici (T, P, umidità)
2. Calcola densità vapor d'acqua
3. Applica modello ITU-R per attenuazione gas
4. Somma contributi di ossigeno e vapor d'acqua

### C. Equazione Radar Completa
1. Calcola distanza 3D radar-target
2. Applica equazione radar base
3. Sottrai attenuazione atmosferica (dB)
4. Confronta con sensibilità minima radar

### D. Rendering 2D
1. Proietta coordinate geografiche su piano 2D
2. Usa proiezione equirettangolare semplificata
3. Scala distanze per visualizzazione
4. Disegna cerchi di copertura e linee LOS

## 3. Implementazione in GolangHo creato un sistema completo di simulazione radar in Go che implementa tutte le caratteristiche richieste. Ecco un riepilogo delle funzionalità:

## Caratteristiche Implementate

### 1. **Geodesia e Geometria Spaziale**
- Calcolo distanze con formula di Haversine
- Conversione coordinate geografiche (lat/lon) in distanze metriche
- Calcolo azimuth e angolo di elevazione
- Supporto per elevazione sul livello del mare

### 2. **Fisica Radar**
- **Equazione radar completa** con attenuazione atmosferica
- Calcolo **orizzonte radar** con fattore di rifrazione K=4/3
- Verifica **Line of Sight** considerando curvatura terrestre
- Calcolo **SNR** (Signal-to-Noise Ratio)

### 3. **Modello Atmosferico**
- Attenuazione per ossigeno e vapor d'acqua (modello ITU-R)
- Dipendenza da temperatura, umidità e pressione
- Calcolo densità vapor d'acqua con formula di Magnus

### 4. **Simulazione Multi-Radar/Multi-Target**
- Supporto per multiple postazioni radar
- Tracking simultaneo di più velivoli
- Risultati dettagliati per ogni coppia radar-target

### 5. **Visualizzazione 2D**
- Proiezione equirettangolare su mappa ASCII
- Rappresentazione posizioni radar e target
- Scala automatica basata su coordinate

## Come Usare

Per eseguire il programma:

```bash
# Salva il codice in un file main.go
go run main.go
```

Il programma simula uno scenario con:
- 2 radar (Roma e Napoli)
- 2 target aerei a diverse altitudini
- Condizioni atmosferiche standard

## Output

Il sistema fornisce:
1. **Parametri atmosferici** utilizzati
2. **Caratteristiche radar** (frequenza, potenza, orizzonte)
3. **Analisi detection** per ogni target:
   - Distanza 3D
   - Azimuth ed elevazione
   - Stato Line of Sight
   - Attenuazione atmosferica
   - Potenza ricevuta e SNR
   - Stato di rilevamento
4. **Mappa tattica ASCII** 2D
5. **Analisi geometrica** distanze inter-radar

## Estensioni Possibili

Il sistema può essere facilmente esteso per:
- Export dati in JSON/CSV
- Rendering grafico con librerie 2D/3D
- Simulazione temporale con movimenti target
- Integrazione con mappe reali (GeoJSON)
- Modelli atmosferici più complessi
- Pattern antenna direzionali
- Jamming ed ECM (Electronic Countermeasures)

Tutte le equazioni sono implementate fedelmente ai modelli fisici standard utilizzati in sistemi radar reali!