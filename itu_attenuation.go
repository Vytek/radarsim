package main

import (
	"math"
)

// P676 implementa ITU-R P.676 per attenuazione da gas atmosferici
type P676 struct {
	Frequency   float64 // GHz
	Pressure    float64 // hPa
	Temperature float64 // Kelvin
	WaterVapor  float64 // g/m³
}

// SpecificAttenuation calcola attenuazione specifica da gas (dB/km)
func (p *P676) SpecificAttenuation() (oxygen, water, total float64) {
	f := p.Frequency
	P := p.Pressure
	T := p.Temperature
	rho := p.WaterVapor

	// Temperatura ridotta
	theta := 300.0 / T

	// Pressione parziale vapor d'acqua (hPa)
	e := rho * T / 216.7

	// Attenuazione da ossigeno (modello semplificato)
	// Risonanze principali O2: 60 GHz e 118.75 GHz
	rp := P * theta

	// Termine non risonante O2
	N_O2 := 6.09 / (rp * math.Pow(theta, 2)) *
		(1 + math.Pow(f/57.0, 2)) *
		math.Pow(f*rp, 2) * 1e-3

	// Risonanza 60 GHz O2 (semplificata)
	f60 := []float64{50.474214, 50.987025, 51.503360, 52.021429,
		52.542418, 53.066934, 53.595775, 54.130025,
		54.671180, 55.221384, 55.783815, 56.264774,
		56.363399, 57.612486, 58.323877, 58.446588,
		59.164204, 59.590983, 60.306056, 60.434778,
		61.150562, 61.800158, 62.411220, 62.486253,
		62.997984, 63.568526, 64.127775, 64.678910,
		65.224078, 65.764779, 66.302096, 66.836834,
		67.369601, 67.900868, 68.431006, 68.960312}

	S_O2 := 0.0
	for _, fi := range f60 {
		delta := 0.0 // larghezza linea semplificata
		if fi < 60 {
			delta = 0.5 * rp * math.Pow(theta, 0.8)
		} else {
			delta = 0.8 * rp * math.Pow(theta, 0.8)
		}
		S_O2 += delta / ((math.Pow(f-fi, 2) + math.Pow(delta, 2)) *
			(math.Pow(f+fi, 2) + math.Pow(delta, 2)))
	}

	oxygen = (7.34e-3*rp*math.Pow(theta, 2)*S_O2 + N_O2) *
		math.Pow(f, 2) * 1e-3

	// Attenuazione da vapor d'acqua (semplificata)
	// Risonanza principale: 22.235 GHz
	f22 := 22.235
	gamma22 := 0.1 * e * math.Pow(theta, 2.5)

	S_H2O := gamma22 / (math.Pow(f-f22, 2) + math.Pow(gamma22, 2))

	// Risonanza 183.31 GHz
	f183 := 183.31
	gamma183 := 0.1 * e * math.Pow(theta, 2.5)
	S_H2O += gamma183 / (math.Pow(f-f183, 2) + math.Pow(gamma183, 2))

	// Contributo continuo
	water = (0.05 + 0.0021*e + S_H2O*3.6) * e * math.Pow(f, 2) *
		math.Pow(theta, 2.5) * 1e-4

	total = oxygen + water
	return
}

// P838 implementa ITU-R P.838 per attenuazione da pioggia
type P838 struct {
	Frequency    float64 // GHz
	Polarization string  // "H" (orizzontale) o "V" (verticale)
	Elevation    float64 // angolo di elevazione (gradi)
	RainRate     float64 // intensità pioggia (mm/h)
}

// Coefficienti k e alpha per P.838
func (p *P838) getCoefficients() (k, alpha float64) {
	f := p.Frequency

	// Coefficienti per polarizzazione orizzontale
	kH := [][]float64{
		{1, -5.33980, -0.10008, 1.13098},
		{2, -0.35351, 1.26970, 0.45400},
		{3, -0.23789, 0.86036, 0.15354},
		{4, -0.94158, 0.64552, 0.16817},
	}

	aH := [][]float64{
		{1, -0.14318, 0.29591, 0.32177},
		{2, -0.29545, 0.77136, 0.67849},
		{3, -0.18363, 0.63297, 0.43386},
		{4, -0.36676, 0.54639, 0.26074},
	}

	// Coefficienti per polarizzazione verticale
	kV := [][]float64{
		{1, -3.80595, 0.56934, 0.81061},
		{2, -3.44965, -0.22911, 0.51059},
		{3, -0.39902, 0.73042, 0.11899},
		{4, 0.50167, 0.07216, 0.27195},
	}

	aV := [][]float64{
		{1, -0.07771, 0.29971, 0.31255},
		{2, 0.56727, -0.20238, 0.33449},
		{3, -0.20238, 0.39570, 0.19408},
		{4, 0.20545, 0.15065, 0.16817},
	}

	logF := math.Log10(f)

	var kCoeff, aCoeff [][]float64
	if p.Polarization == "V" {
		kCoeff = kV
		aCoeff = aV
	} else {
		kCoeff = kH
		aCoeff = aH
	}

	// Calcolo k
	logK := 0.0
	for _, row := range kCoeff {
		logK += row[1] * math.Exp(-math.Pow((logF-row[2])/row[3], 2))
	}
	k = math.Pow(10, logK)

	// Calcolo alpha
	alpha = 0.0
	for _, row := range aCoeff {
		alpha += row[1] * math.Exp(-math.Pow((logF-row[2])/row[3], 2))
	}

	return
}

// SpecificAttenuation calcola attenuazione specifica da pioggia (dB/km)
func (p *P838) SpecificAttenuation() float64 {
	k, alpha := p.getCoefficients()
	gammaR := k * math.Pow(p.RainRate, alpha)
	return gammaR
}

// PathAttenuation calcola attenuazione su percorso inclinato
func (p *P838) PathAttenuation(pathLength float64) float64 {
	gammaR := p.SpecificAttenuation()

	// Fattore di riduzione per percorso inclinato
	theta := p.Elevation
	var r float64
	if theta >= 5 {
		r = 1.0 / (1 + pathLength*gammaR/36000)
	} else {
		// Per angoli bassi, fattore ridotto
		r = 1.0 / (1 + pathLength*gammaR/54000)
	}

	return gammaR * pathLength * r
}

// P452 struttura semplificata per path loss
type P452 struct {
	Frequency   float64   // MHz
	Distance    float64   // km
	TxHeight    float64   // m
	RxHeight    float64   // m
	PathProfile []float64 // altezze terreno (m)
	TimePercent float64   // percentuale tempo (%)
}

// BasicPathLoss calcola path loss base (free space + diffrazione)
func (p *P452) BasicPathLoss() float64 {
	f := p.Frequency
	d := p.Distance

	// Free space path loss
	Lbf := 32.45 + 20*math.Log10(f) + 20*math.Log10(d)

	// Diffrazione (modello semplificato knife-edge)
	// Trova ostacolo più alto
	maxHeight := 0.0
	if len(p.PathProfile) > 0 {
		for _, h := range p.PathProfile {
			if h > maxHeight {
				maxHeight = h
			}
		}
	}

	// Calcola parametro di diffrazione v
	h := maxHeight - (p.TxHeight+p.RxHeight)/2
	if h > 0 {
		lambda := 300.0 / f // wavelength in m
		v := h * math.Sqrt(2*p.Distance*1000/(lambda*d*1000))

		// Attenuazione diffrazione (dB)
		var Ld float64
		if v > -0.78 {
			Ld = 6.9 + 20*math.Log10(math.Sqrt(math.Pow(v-0.1, 2)+1)+v-0.1)
		} else {
			Ld = 0
		}

		return Lbf + Ld
	}

	return Lbf
}

// TroposcatterLoss calcola perdita troposferica
func (p *P452) TroposcatterLoss() float64 {
	f := p.Frequency
	d := p.Distance

	// Modello semplificato
	Lbs := 190 + 20*math.Log10(d) + 30*math.Log10(f) - 0.3*p.TimePercent

	return Lbs
}

// Esempio di utilizzo
/*
func main() {
	fmt.Println("=== ITU-R P.676: Attenuazione da Gas ===")
	p676 := P676{
		Frequency:   30.0,  // 30 GHz
		Pressure:    1013.25, // livello mare
		Temperature: 288.15,  // 15°C
		WaterVapor:  7.5,     // g/m³
	}

	o2, h2o, total := p676.SpecificAttenuation()
	fmt.Printf("Frequenza: %.2f GHz\n", p676.Frequency)
	fmt.Printf("Attenuazione O₂: %.4f dB/km\n", o2)
	fmt.Printf("Attenuazione H₂O: %.4f dB/km\n", h2o)
	fmt.Printf("Attenuazione totale: %.4f dB/km\n\n", total)

	fmt.Println("=== ITU-R P.838: Attenuazione da Pioggia ===")
	p838 := P838{
		Frequency:    12.0, // 12 GHz (Ku band)
		Polarization: "H",
		Elevation:    30.0, // gradi
		RainRate:     50.0, // mm/h (pioggia forte)
	}

	gammaR := p838.SpecificAttenuation()
	pathLen := 10.0 // km
	pathAtt := p838.PathAttenuation(pathLen)

	fmt.Printf("Frequenza: %.2f GHz\n", p838.Frequency)
	fmt.Printf("Intensità pioggia: %.1f mm/h\n", p838.RainRate)
	fmt.Printf("Attenuazione specifica: %.4f dB/km\n", gammaR)
	fmt.Printf("Attenuazione su %g km: %.2f dB\n\n", pathLen, pathAtt)

	fmt.Println("=== ITU-R P.452: Path Loss ===")
	p452 := P452{
		Frequency:   2000.0, // 2 GHz
		Distance:    50.0,   // km
		TxHeight:    30.0,   // m
		RxHeight:    10.0,   // m
		PathProfile: []float64{0, 50, 100, 80, 40, 0}, // profilo terreno
		TimePercent: 50.0,   // 50% del tempo
	}

	basicLoss := p452.BasicPathLoss()
	tropoLoss := p452.TroposcatterLoss()

	fmt.Printf("Frequenza: %.0f MHz\n", p452.Frequency)
	fmt.Printf("Distanza: %.1f km\n", p452.Distance)
	fmt.Printf("Path loss base (free space + diffrazione): %.2f dB\n", basicLoss)
	fmt.Printf("Path loss troposferico: %.2f dB\n", tropoLoss)

	// Esempio combinato
	fmt.Println("\n=== Esempio Link Budget Completo ===")
	freq := 30.0 // GHz
	dist := 10.0 // km

	// Gas
	p676_ex := P676{freq, 1013.25, 288.15, 7.5}
	_, _, gasAtt := p676_ex.SpecificAttenuation()
	totalGasAtt := gasAtt * dist

	// Pioggia
	p838_ex := P838{freq, "V", 45.0, 25.0}
	rainAtt := p838_ex.PathAttenuation(dist)

	// Free space
	fspl := 32.45 + 20*math.Log10(freq*1000) + 20*math.Log10(dist)

	totalLoss := fspl + totalGasAtt + rainAtt

	fmt.Printf("Link a %.0f GHz, distanza %.0f km:\n", freq, dist)
	fmt.Printf("- Free Space Path Loss: %.2f dB\n", fspl)
	fmt.Printf("- Attenuazione gas: %.2f dB\n", totalGasAtt)
	fmt.Printf("- Attenuazione pioggia: %.2f dB\n", rainAtt)
	fmt.Printf("- PERDITA TOTALE: %.2f dB\n", totalLoss)
}
*/
