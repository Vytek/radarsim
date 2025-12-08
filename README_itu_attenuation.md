# Raccomandazioni ITU-R per Attenuazione Atmosferica

Le raccomandazioni ITU-R forniscono modelli standardizzati per calcolare l'attenuazione dei segnali radio causata dall'atmosfera. Ecco le principali:

## Raccomandazioni Principali

**ITU-R P.676** - Attenuazione da gas atmosferici (ossigeno e vapor d'acqua)

- Calcola l'attenuazione specifica γ (dB/km) per frequenze fino a 1000 GHz
- Considera ossigeno (O₂) e vapor d'acqua (H₂O)
- Dipende da frequenza, pressione, temperatura e densità di vapor d'acqua
- Essenziale per collegamenti sopra 10 GHz

**ITU-R P.838** - Attenuazione da pioggia

- Modello per calcolare l'attenuazione specifica causata dalla pioggia
- Usa coefficienti k e α dipendenti dalla frequenza e polarizzazione
- Formula: γ_R = k × R^α (dove R è l'intensità di pioggia in mm/h)
- Valido da 1 a 1000 GHz

**ITU-R P.452** - Predizione di path loss

- Modello completo per collegamenti punto-punto terrestri
- Include diffrazione, troposfera, effetti del terreno
- Considera statistiche climatiche a lungo termine
- Più complesso, integra vari fenomeni di propagazione

## Implementazione in Go

Vedi file itu_attenuation.go

## Note sull'implementazione

**Semplificazioni adottate:**

- P.676: implementazione delle risonanze principali O₂ (60 GHz) e H₂O (22.235 GHz), formula completa include molte più linee spettrali
- P.838: usa fitting polynomial dei coefficienti, la raccomandazione ITU fornisce tabelle dettagliate
- P.452: versione molto semplificata, il modello completo include troposfera, ducting, diffrazione multipla

**Per implementazione production:**

- Usa le tabelle complete ITU-R ufficiali
- Implementa tutti i casi limite e validazioni
- Considera variazioni stagionali e climatiche
- Include modelli statistici per disponibilità del link

**Applicazioni pratiche:**

- Design collegamenti satellitari (soprattutto Ka/V band)
- Reti 5G mmWave (24-40 GHz)
- Collegamenti punto-punto terrestri
- Link budget per sistemi di comunicazione critica