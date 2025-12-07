package main

import (
	"fmt"
	"math"
	"strings"
)

// Costanti fisiche
const (
	EarthRadius      = 6371000.0   // metri
	SpeedOfLight     = 299792458.0 // m/s
	RefractionFactor = 4.0 / 3.0   // Fattore K per rifrazione atmosferica
	AtmPressure      = 1013.25     // hPa (livello del mare)
)

// Coordinate geografiche
type GeoCoord struct {
	Latitude  float64 // gradi
	Longitude float64 // gradi
	Elevation float64 // metri sul livello del mare
}

// Parametri atmosferici
type AtmosphericConditions struct {
	Temperature float64 // Celsius
	Humidity    float64 // % (0-100)
	Pressure    float64 // hPa
}

// Configurazione radar
type RadarStation struct {
	ID             string
	Position       GeoCoord
	Frequency      float64 // Hz
	TxPower        float64 // Watt
	AntennaGain    float64 // dBi
	MinDetectPower float64 // Watt (sensibilità minima)
	MaxRange       float64 // metri
}

// Bersaglio aereo
type AircraftTarget struct {
	ID       string
	Position GeoCoord
	RCS      float64 // m² (Radar Cross Section)
	Velocity float64 // m/s
	Heading  float64 // gradi
}

// Risultato detection
type DetectionResult struct {
	RadarID       string
	TargetID      string
	Distance      float64
	Azimuth       float64
	Elevation     float64
	ReceivedPower float64
	SNR           float64
	IsDetected    bool
	IsLOS         bool
	Attenuation   float64
}

// ========== FUNZIONI GEODETICHE ==========

// Converte gradi in radianti
func deg2rad(deg float64) float64 {
	return deg * math.Pi / 180.0
}

// Converte radianti in gradi
func rad2deg(rad float64) float64 {
	return rad * 180.0 / math.Pi
}

// Calcola distanza orizzontale usando formula Haversine
func haversineDistance(coord1, coord2 GeoCoord) float64 {
	lat1, lon1 := deg2rad(coord1.Latitude), deg2rad(coord1.Longitude)
	lat2, lon2 := deg2rad(coord2.Latitude), deg2rad(coord2.Longitude)

	dlat := lat2 - lat1
	dlon := lon2 - lon1

	a := math.Sin(dlat/2)*math.Sin(dlat/2) +
		math.Cos(lat1)*math.Cos(lat2)*
			math.Sin(dlon/2)*math.Sin(dlon/2)

	c := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))

	return EarthRadius * c
}

// Calcola distanza orizzontale usando formula di Vincenty (più accurata)
// Utilizza il modello ellissoidale WGS-84 della Terra
func vincentyDistance(coord1, coord2 GeoCoord) float64 {
	// Costanti ellissoide WGS-84
	const (
		a = 6378137.0         // semiasse maggiore (metri)
		b = 6356752.314245    // semiasse minore (metri)
		f = 1 / 298.257223563 // appiattimento
	)

	lat1, lon1 := deg2rad(coord1.Latitude), deg2rad(coord1.Longitude)
	lat2, lon2 := deg2rad(coord2.Latitude), deg2rad(coord2.Longitude)

	L := lon2 - lon1
	U1 := math.Atan((1 - f) * math.Tan(lat1))
	U2 := math.Atan((1 - f) * math.Tan(lat2))

	sinU1 := math.Sin(U1)
	cosU1 := math.Cos(U1)
	sinU2 := math.Sin(U2)
	cosU2 := math.Cos(U2)

	lambda := L
	lambdaP := 2 * math.Pi
	iterLimit := 100
	var cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma float64

	for math.Abs(lambda-lambdaP) > 1e-12 && iterLimit > 0 {
		sinLambda := math.Sin(lambda)
		cosLambda := math.Cos(lambda)

		sinSigma = math.Sqrt((cosU2*sinLambda)*(cosU2*sinLambda) +
			(cosU1*sinU2-sinU1*cosU2*cosLambda)*(cosU1*sinU2-sinU1*cosU2*cosLambda))

		if sinSigma == 0 {
			return 0 // Punti coincidenti
		}

		cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
		sigma = math.Atan2(sinSigma, cosSigma)
		sinAlpha := cosU1 * cosU2 * sinLambda / sinSigma
		cosSqAlpha = 1 - sinAlpha*sinAlpha

		cos2SigmaM = 0.0
		if cosSqAlpha != 0 {
			cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha
		}

		C := f / 16 * cosSqAlpha * (4 + f*(4-3*cosSqAlpha))
		lambdaP = lambda
		lambda = L + (1-C)*f*sinAlpha*
			(sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))

		iterLimit--
	}

	if iterLimit == 0 {
		// Fallback su Haversine se non converge
		return haversineDistance(coord1, coord2)
	}

	uSq := cosSqAlpha * (a*a - b*b) / (b * b)
	A := 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	B := uSq / 1024 * (256 + uSq*(-128+uSq*(74-47*uSq)))

	deltaSigma := B * sinSigma * (cos2SigmaM + B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
		B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))

	distance := b * A * (sigma - deltaSigma)

	return distance
}

// Calcola distanza 3D includendo differenza di elevazione
func distance3D(coord1, coord2 GeoCoord) float64 {
	hDist := haversineDistance(coord1, coord2)
	vDist := coord2.Elevation - coord1.Elevation
	return math.Sqrt(hDist*hDist + vDist*vDist)
}

// Calcola distanza 3D usando Vincenty (più accurata)
func distance3DVincenty(coord1, coord2 GeoCoord) float64 {
	hDist := vincentyDistance(coord1, coord2)
	vDist := coord2.Elevation - coord1.Elevation
	return math.Sqrt(hDist*hDist + vDist*vDist)
}

// Calcola azimuth (bearing) da coord1 a coord2
func calculateAzimuth(coord1, coord2 GeoCoord) float64 {
	lat1, lon1 := deg2rad(coord1.Latitude), deg2rad(coord1.Longitude)
	lat2, lon2 := deg2rad(coord2.Latitude), deg2rad(coord2.Longitude)

	dlon := lon2 - lon1

	y := math.Sin(dlon) * math.Cos(lat2)
	x := math.Cos(lat1)*math.Sin(lat2) -
		math.Sin(lat1)*math.Cos(lat2)*math.Cos(dlon)

	azimuth := math.Atan2(y, x)
	return math.Mod(rad2deg(azimuth)+360, 360)
}

// Calcola angolo di elevazione
func calculateElevation(coord1, coord2 GeoCoord) float64 {
	hDist := haversineDistance(coord1, coord2)
	vDist := coord2.Elevation - coord1.Elevation

	if hDist == 0 {
		if vDist > 0 {
			return 90.0
		} else if vDist < 0 {
			return -90.0
		}
		return 0.0
	}

	return rad2deg(math.Atan2(vDist, hDist))
}

// ========== ORIZZONTE RADAR ==========

// Calcola orizzonte radar per una data altezza
func radarHorizon(heightMeters float64) float64 {
	return math.Sqrt(2.0 * RefractionFactor * EarthRadius * heightMeters)
}

// Verifica Line of Sight considerando curvatura terrestre
func hasLineOfSight(radar, target GeoCoord) bool {
	horizonRadar := radarHorizon(radar.Elevation)
	horizonTarget := radarHorizon(target.Elevation)

	hDist := haversineDistance(radar, target)

	// Distanza massima LOS è la somma dei due orizzonti
	maxLOSDistance := horizonRadar + horizonTarget

	return hDist <= maxLOSDistance
}

// ========== ATTENUAZIONE ATMOSFERICA ==========

// Calcola densità vapor d'acqua (g/m³)
func waterVaporDensity(temp, humidity float64) float64 {
	// Pressione di saturazione del vapor d'acqua (Magnus formula)
	es := 6.1078 * math.Exp((17.27*temp)/(temp+237.3))

	// Pressione parziale vapor d'acqua
	e := (humidity / 100.0) * es

	// Densità vapor d'acqua
	return (216.7 * e) / (temp + 273.15)
}

// Calcola attenuazione atmosferica (dB/km) - modello semplificato ITU-R
func atmosphericAttenuation(freqHz float64, conditions AtmosphericConditions, distance float64) float64 {
	freqGHz := freqHz / 1e9
	temp := conditions.Temperature
	humidity := conditions.Humidity

	// Densità vapor d'acqua
	rho := waterVaporDensity(temp, humidity)

	// Attenuazione per ossigeno (approssimazione)
	gammaO2 := (7.19e-3 * freqGHz * freqGHz * conditions.Pressure) /
		((freqGHz*freqGHz + 0.34) * (freqGHz*freqGHz + 0.27*conditions.Pressure*conditions.Pressure))

	// Attenuazione per vapor d'acqua (approssimazione)
	term1 := (freqGHz - 22.2) * (freqGHz - 22.2)
	term2 := (freqGHz - 183.3) * (freqGHz - 183.3)
	gammaH2O := (0.05 + 0.0021*rho + 3.6e-4*rho*rho) * freqGHz * freqGHz / (term1 + 8.5 + term2 + 9.0)

	// Attenuazione totale per unità di lunghezza
	gammaTotal := gammaO2 + gammaH2O

	// Attenuazione totale sul percorso (dB)
	return gammaTotal * (distance / 1000.0) // converti distanza in km
}

// ========== EQUAZIONE RADAR ==========

// Calcola potenza ricevuta usando equazione radar
func radarEquation(radar RadarStation, target AircraftTarget, distance float64, attenuation float64) float64 {
	// Lunghezza d'onda
	wavelength := SpeedOfLight / radar.Frequency

	// Conversione guadagno da dBi a lineare
	gainLinear := math.Pow(10, radar.AntennaGain/10.0)

	// Equazione radar (forma semplificata, Gt = Gr)
	numerator := radar.TxPower * gainLinear * gainLinear * wavelength * wavelength * target.RCS
	denominator := math.Pow(4*math.Pi, 3) * math.Pow(distance, 4)

	receivedPower := numerator / denominator

	// Applica attenuazione atmosferica (converti dB in fattore lineare)
	attenuationFactor := math.Pow(10, -attenuation/10.0)
	receivedPower *= attenuationFactor

	return receivedPower
}

// Calcola SNR (Signal to Noise Ratio)
func calculateSNR(receivedPower, noisePower float64) float64 {
	if noisePower <= 0 {
		noisePower = 1e-15 // Valore minimo per evitare divisione per zero
	}
	return 10 * math.Log10(receivedPower/noisePower)
}

// ========== SIMULAZIONE ==========

// Simula detection di un target da un radar
func simulateDetection(radar RadarStation, target AircraftTarget, conditions AtmosphericConditions) DetectionResult {
	result := DetectionResult{
		RadarID:  radar.ID,
		TargetID: target.ID,
	}

	// Calcola geometria
	result.Distance = distance3D(radar.Position, target.Position)
	result.Azimuth = calculateAzimuth(radar.Position, target.Position)
	result.Elevation = calculateElevation(radar.Position, target.Position)

	// Verifica Line of Sight
	result.IsLOS = hasLineOfSight(radar.Position, target.Position)

	if !result.IsLOS {
		result.IsDetected = false
		return result
	}

	// Verifica range massimo
	if result.Distance > radar.MaxRange {
		result.IsDetected = false
		return result
	}

	// Calcola attenuazione atmosferica
	result.Attenuation = atmosphericAttenuation(radar.Frequency, conditions, result.Distance)

	// Calcola potenza ricevuta
	result.ReceivedPower = radarEquation(radar, target, result.Distance, result.Attenuation)

	// Calcola SNR
	result.SNR = calculateSNR(result.ReceivedPower, radar.MinDetectPower)

	// Detection avviene se potenza ricevuta > sensibilità minima
	result.IsDetected = result.ReceivedPower >= radar.MinDetectPower

	return result
}

// ========== VISUALIZZAZIONE 2D ==========

// Proiezione coordinate per visualizzazione 2D
type Point2D struct {
	X float64
	Y float64
}

// Converte coordinate geografiche in punto 2D (proiezione equirettangolare)
func geoToPoint2D(coord GeoCoord, center GeoCoord, scale float64) Point2D {
	// Offset dal centro
	dLat := coord.Latitude - center.Latitude
	dLon := coord.Longitude - center.Longitude

	// Proiezione equirettangolare con correzione latitudine
	cosLat := math.Cos(deg2rad(center.Latitude))

	return Point2D{
		X: dLon * cosLat * scale,
		Y: dLat * scale,
	}
}

// Genera rappresentazione ASCII della situazione tattica
func generateTacticalMap(radars []RadarStation, targets []AircraftTarget, center GeoCoord, width, height int) string {
	// Crea griglia
	grid := make([][]rune, height)
	for i := range grid {
		grid[i] = make([]rune, width)
		for j := range grid[i] {
			grid[i][j] = '.'
		}
	}

	// Scala per conversione coordinate
	scale := float64(width) / 4.0

	// Disegna radar
	for _, radar := range radars {
		pt := geoToPoint2D(radar.Position, center, scale)
		x := int(pt.X) + width/2
		y := height/2 - int(pt.Y)

		if x >= 0 && x < width && y >= 0 && y < height {
			grid[y][x] = 'R'
		}
	}

	// Disegna target
	for _, target := range targets {
		pt := geoToPoint2D(target.Position, center, scale)
		x := int(pt.X) + width/2
		y := height/2 - int(pt.Y)

		if x >= 0 && x < width && y >= 0 && y < height {
			grid[y][x] = 'T'
		}
	}

	// Converti griglia in stringa
	var builder strings.Builder
	for _, row := range grid {
		builder.WriteString(string(row))
		builder.WriteString("\n")
	}

	return builder.String()
}

// ========== FUNZIONE PRINCIPALE ==========

func main() {
	fmt.Println("=== SISTEMA DI SIMULAZIONE RADAR ===\n")

	// Condizioni atmosferiche
	conditions := AtmosphericConditions{
		Temperature: 15.0,    // °C
		Humidity:    60.0,    // %
		Pressure:    1013.25, // hPa
	}

	// Configura radar (esempio: Roma)
	radar1 := RadarStation{
		ID: "RADAR-01",
		Position: GeoCoord{
			Latitude:  41.9028,
			Longitude: 12.4964,
			Elevation: 100.0, // 100m slm
		},
		Frequency:      3e9,     // 3 GHz (banda S)
		TxPower:        1000000, // 1 MW
		AntennaGain:    30.0,    // 30 dBi
		MinDetectPower: 1e-12,   // -120 dBm
		MaxRange:       400000,  // 400 km
	}

	// Configura secondo radar (esempio: Napoli)
	radar2 := RadarStation{
		ID: "RADAR-02",
		Position: GeoCoord{
			Latitude:  40.8518,
			Longitude: 14.2681,
			Elevation: 150.0,
		},
		Frequency:      3e9,
		TxPower:        1000000,
		AntennaGain:    30.0,
		MinDetectPower: 1e-12,
		MaxRange:       400000,
	}

	radars := []RadarStation{radar1, radar2}

	// Configura target aerei
	target1 := AircraftTarget{
		ID: "TGT-001",
		Position: GeoCoord{
			Latitude:  42.0,
			Longitude: 12.8,
			Elevation: 10000.0, // 10km altitudine
		},
		RCS:      10.0,  // 10 m² (caccia tattico)
		Velocity: 250.0, // 900 km/h
		Heading:  180.0,
	}

	target2 := AircraftTarget{
		ID: "TGT-002",
		Position: GeoCoord{
			Latitude:  41.5,
			Longitude: 13.5,
			Elevation: 8000.0,
		},
		RCS:      50.0, // 50 m² (aereo commerciale)
		Velocity: 230.0,
		Heading:  270.0,
	}

	targets := []AircraftTarget{target1, target2}

	// Esegui simulazione
	fmt.Println("CONDIZIONI ATMOSFERICHE:")
	fmt.Printf("  Temperatura: %.1f°C\n", conditions.Temperature)
	fmt.Printf("  Umidità: %.1f%%\n", conditions.Humidity)
	fmt.Printf("  Pressione: %.2f hPa\n\n", conditions.Pressure)

	// Simula detection per ogni combinazione radar-target
	for _, radar := range radars {
		fmt.Printf("RADAR: %s (Lat: %.4f°, Lon: %.4f°, Alt: %.0fm)\n",
			radar.ID, radar.Position.Latitude, radar.Position.Longitude, radar.Position.Elevation)
		fmt.Printf("  Frequenza: %.2f GHz\n", radar.Frequency/1e9)
		fmt.Printf("  Potenza TX: %.0f kW\n", radar.TxPower/1000)
		fmt.Printf("  Orizzonte Radar: %.1f km\n\n", radarHorizon(radar.Position.Elevation)/1000)

		for _, target := range targets {
			result := simulateDetection(radar, target, conditions)

			fmt.Printf("  TARGET: %s\n", result.TargetID)
			fmt.Printf("    Posizione: Lat %.4f°, Lon %.4f°, Alt %.0fm\n",
				target.Position.Latitude, target.Position.Longitude, target.Position.Elevation)
			fmt.Printf("    RCS: %.1f m²\n", target.RCS)
			fmt.Printf("    Distanza: %.2f km\n", result.Distance/1000)
			fmt.Printf("    Azimuth: %.2f°\n", result.Azimuth)
			fmt.Printf("    Elevazione: %.2f°\n", result.Elevation)
			fmt.Printf("    Line of Sight: %v\n", result.IsLOS)
			fmt.Printf("    Attenuazione atmosferica: %.2f dB\n", result.Attenuation)
			fmt.Printf("    Potenza ricevuta: %.2e W (%.2f dBm)\n",
				result.ReceivedPower, 10*math.Log10(result.ReceivedPower/1e-3))
			fmt.Printf("    SNR: %.2f dB\n", result.SNR)
			fmt.Printf("    RILEVATO: %v\n\n", result.IsDetected)
		}
	}

	// Genera mappa tattica
	center := GeoCoord{
		Latitude:  41.5,
		Longitude: 13.0,
		Elevation: 0,
	}

	fmt.Println("MAPPA TATTICA (R=Radar, T=Target):")
	fmt.Println(generateTacticalMap(radars, targets, center, 60, 30))

	// Calcoli aggiuntivi
	fmt.Println("ANALISI GEOMETRICA:")
	for i := 0; i < len(radars); i++ {
		for j := i + 1; j < len(radars); j++ {
			dist := haversineDistance(radars[i].Position, radars[j].Position)
			fmt.Printf("  Distanza %s - %s: %.2f km\n",
				radars[i].ID, radars[j].ID, dist/1000)
		}
	}
}
