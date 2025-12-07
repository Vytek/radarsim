// radar_sim_go.go
// Implementazione di base in Go per una simulazione radar che
// considera lat/lon/elev, orizzonte geometrico/radar (k-factor),
// equazione radar monostatica con attenuazione atmosferica semplice,
// e generazione di poligoni di copertura GeoJSON (2D vettoriale).
//
// Questo file è pensato come punto di partenza: include strutture,
// funzioni geodetiche (WGS84), equazione radar, algoritmo per
// calcolare r(azimuth) dove SNR supera una soglia, e esportazione GeoJSON.
// Miglioramenti possibili: integrazione DEM, modelli ITU per attenuazione,
// ray-tracing, pattern d'antenna, processamento impulsi, multi-threading.

package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

const VERSION = "0.0.2"

// ------------------ Costanti WGS84 ------------------
const (
	WGS84_A  = 6378137.0               // semi-major axis (m)
	WGS84_F  = 1.0 / 298.257223563     // flattening
	WGS84_E2 = WGS84_F * (2 - WGS84_F) // eccentricity^2
	C        = 299792458.0             // speed of light (m/s)
	K_B      = 1.380649e-23            // Boltzmann constant
)

// ------------------ Tipi principali ------------------

// Radar rappresenta una postazione radar
type Radar struct {
	ID           string
	Lat          float64 // degrees
	Lon          float64 // degrees
	Height       float64 // meters above MSL
	Pt           float64 // transmitted power (W)
	Gain         float64 // linear
	Freq         float64 // Hz
	Sigma        float64 // assumed RCS (m^2) for coverage calc
	Loss         float64 // linear system losses
	Ts           float64 // system noise temperature (K)
	Bandwidth    float64 // receiver bandwidth (Hz)
	SNRThreshold float64 // linear SNR threshold
	KFactor      float64 // effective earth radius factor (e.g., 4/3)
}

// Target rappresenta un bersaglio (usato qui per test)
type Target struct {
	Lat    float64
	Lon    float64
	Height float64
	RCS    float64
}

// Environment parametri atmosferici semplici
type Environment struct {
	TempC       float64 // Celsius
	RH          float64 // relative humidity (0-100)
	PressureHPa float64
}

// GeoJSON structures (minimal)
type GeoJSONFeatureCollection struct {
	Type     string           `json:"type"`
	Features []GeoJSONFeature `json:"features"`
}

type GeoJSONFeature struct {
	Type       string                 `json:"type"`
	Properties map[string]interface{} `json:"properties"`
	Geometry   GeoJSONGeometry        `json:"geometry"`
}

type GeoJSONGeometry struct {
	Type        string        `json:"type"`
	Coordinates [][][]float64 `json:"coordinates"` // only for Polygon
}

// ------------------ Geodetica: lat/lon <-> ECEF / ENU ------------------

func deg2rad(d float64) float64 { return d * math.Pi / 180.0 }
func rad2deg(r float64) float64 { return r * 180.0 / math.Pi }

// LatLonAlt to ECEF (meters), ECEF = Earth-Centered, Earth-Fixed
func llhToEcef(latDeg, lonDeg, h float64) (x, y, z float64) {
	phi := deg2rad(latDeg)
	lambda := deg2rad(lonDeg)
	N := WGS84_A / math.Sqrt(1-WGS84_E2*math.Sin(phi)*math.Sin(phi))
	x = (N + h) * math.Cos(phi) * math.Cos(lambda)
	y = (N + h) * math.Cos(phi) * math.Sin(lambda)
	z = ((1-WGS84_E2)*N + h) * math.Sin(phi)
	return
}

// ECEF -> vector difference, then rotate to ENU at lat0, lon0; ENU = East, North, Up
func ecefToEnu(dx, dy, dz, lat0Deg, lon0Deg float64) (e, n, u float64) {
	phi := deg2rad(lat0Deg)
	lambda := deg2rad(lon0Deg)
	// rotation matrix components
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	sinlam := math.Sin(lambda)
	coslam := math.Cos(lambda)
	// ENU = R * dX
	e = -sinlam*dx + coslam*dy
	n = -coslam*sinphi*dx - sinphi*sinlam*dy + cosphi*dz
	u = cosphi*coslam*dx + cosphi*sinlam*dy + sinphi*dz
	return
}

// Haversine distance (meters) and bearing (deg) from p1 to p2
func haversineDistance(lat1, lon1, lat2, lon2 float64) (d float64, bearingDeg float64) {
	phi1 := deg2rad(lat1)
	phi2 := deg2rad(lat2)
	dphi := deg2rad(lat2 - lat1)
	dlambda := deg2rad(lon2 - lon1)
	a := math.Sin(dphi/2)*math.Sin(dphi/2) + math.Cos(phi1)*math.Cos(phi2)*math.Sin(dlambda/2)*math.Sin(dlambda/2)
	c := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))
	d = WGS84_A * c
	// initial bearing
	y := math.Sin(dlambda) * math.Cos(phi2)
	x := math.Cos(phi1)*math.Sin(phi2) - math.Sin(phi1)*math.Cos(phi2)*math.Cos(dlambda)
	bearing := math.Atan2(y, x)
	bearingDeg = math.Mod(rad2deg(bearing)+360.0, 360.0)
	return
}

// Destination point given start lat/lon, bearing (deg) and distance (meters)
func destinationPoint(lat1, lon1, bearingDeg, distance float64) (lat2, lon2 float64) {
	phi1 := deg2rad(lat1)
	lambda1 := deg2rad(lon1)
	theta := deg2rad(bearingDeg)
	delta := distance / WGS84_A
	phi2 := math.Asin(math.Sin(phi1)*math.Cos(delta) + math.Cos(phi1)*math.Sin(delta)*math.Cos(theta))
	lambda2 := lambda1 + math.Atan2(math.Sin(theta)*math.Sin(delta)*math.Cos(phi1), math.Cos(delta)-math.Sin(phi1)*math.Sin(phi2))
	lat2 = rad2deg(phi2)
	lon2 = rad2deg(lambda2)
	return
}

// Vincenty inverse formula: returns distance in meters and initial bearing in degrees
func vincentyInverse(lat1, lon1, lat2, lon2 float64) (distance, bearing float64) {
	phi1 := deg2rad(lat1)
	phi2 := deg2rad(lat2)
	L := deg2rad(lon2 - lon1)
	u1 := math.Atan((1 - WGS84_F) * math.Tan(phi1))
	u2 := math.Atan((1 - WGS84_F) * math.Tan(phi2))
	sinU1, cosU1 := math.Sin(u1), math.Cos(u1)
	sinU2, cosU2 := math.Sin(u2), math.Cos(u2)

	lambda := L
	lambdaP := 2 * math.Pi
	iterLimit := 100
	var sinSigma, cosSigma, sigma, sinAlpha, cos2Alpha, cos2SigmaM float64

	for iter := 0; math.Abs(lambda-lambdaP) > 1e-12 && iter < iterLimit; iter++ {
		sinLambda := math.Sin(lambda)
		cosLambda := math.Cos(lambda)
		sinSigma = math.Sqrt(math.Pow(cosU2*sinLambda, 2) + math.Pow(cosU1*sinU2-sinU1*cosU2*cosLambda, 2))
		if sinSigma == 0 {
			return 0, 0
		}
		cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
		sigma = math.Atan2(sinSigma, cosSigma)
		sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
		cos2Alpha = 1 - sinAlpha*sinAlpha
		if cos2Alpha == 0 {
			cos2SigmaM = 0
		} else {
			cos2SigmaM = cosSigma - 2*sinU1*sinU2/cos2Alpha
		}
		C := WGS84_F / 16 * cos2Alpha * (4 + WGS84_F*(4-3*cos2Alpha))
		lambdaP = lambda
		lambda = L + (1-C)*WGS84_F*sinAlpha*(sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
	}

	a := WGS84_A
	b := WGS84_A * (1 - WGS84_F)
	uSq := cos2Alpha * (a*a - b*b) / (b * b)
	A := 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	B := uSq / 1024 * (256 + uSq*(-128+uSq*(74-47*uSq)))
	deltaSigma := B * sinSigma * (cos2SigmaM + B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))
	distance = b * A * (sigma - deltaSigma)

	y := cosU2 * math.Sin(lambda)
	x := cosU1*sinU2 - sinU1*cosU2*math.Cos(lambda)
	initialBearing := math.Mod(rad2deg(math.Atan2(y, x))+360, 360)
	return distance, initialBearing
}

// Vincenty direct formula: returns destination lat/lon given start point, bearing (deg) and distance (m)
func destinationPointVincenty(lat1, lon1, bearingDeg, distance float64) (lat2, lon2 float64) {
	// Use Vincenty direct solution (simplified for small distances)
	alpha1 := deg2rad(bearingDeg)
	s := distance
	a := WGS84_A
	f := WGS84_F
	b := a * (1 - f)
	phi1 := deg2rad(lat1)
	lambda1 := deg2rad(lon1)
	u1 := math.Atan((1 - f) * math.Tan(phi1))
	sinU1, cosU1 := math.Sin(u1), math.Cos(u1)
	sinAlpha1, cosAlpha1 := math.Sin(alpha1), math.Cos(alpha1)
	sigma1 := math.Atan2(math.Tan(u1), cosAlpha1)
	sinAlpha := cosU1 * sinAlpha1
	cos2Alpha := 1 - sinAlpha*sinAlpha
	uSq := cos2Alpha * (a*a - b*b) / (b * b)
	A := 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	B := uSq / 1024 * (256 + uSq*(-128+uSq*(74-47*uSq)))
	sigma := s / (b * A)
	sigmaP := 2 * math.Pi
	iterLimit := 100
	var cos2SigmaM float64
	for iter := 0; math.Abs(sigma-sigmaP) > 1e-12 && iter < iterLimit; iter++ {
		cos2SigmaM = math.Cos(2*sigma1 + sigma)
		deltaSigma := B * math.Sin(sigma) * (cos2SigmaM + B/4*(math.Cos(sigma)*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*math.Sin(sigma)*math.Sin(sigma))*(-3+4*cos2SigmaM*cos2SigmaM)))
		sigmaP = sigma
		sigma = s/(b*A) + deltaSigma
	}
	phi2 := math.Atan2(sinU1*math.Cos(sigma)+cosU1*math.Sin(sigma)*cosAlpha1, (1-f)*math.Sqrt(sinAlpha*sinAlpha+math.Pow(sinU1*math.Sin(sigma)-cosU1*math.Cos(sigma)*cosAlpha1, 2)))
	lambda := math.Atan2(math.Sin(sigma)*sinAlpha1, cosU1*math.Cos(sigma)-sinU1*math.Sin(sigma)*cosAlpha1)
	C := f / 16 * cos2Alpha * (4 + f*(4-3*cos2Alpha))
	L := lambda - (1-C)*f*sinAlpha*(sigma+C*math.Sin(sigma)*(cos2SigmaM+C*math.Cos(sigma)*(-1+2*cos2SigmaM*cos2SigmaM)))
	lon2 = rad2deg(lambda1 + L)
	lat2 = rad2deg(phi2)
	return
}

// Compute 3D distance from radar to a point at given lat/lon/height
func distance3D(lat1, lon1, h1, lat2, lon2, h2 float64) float64 {
	x1, y1, z1 := llhToEcef(lat1, lon1, h1)
	x2, y2, z2 := llhToEcef(lat2, lon2, h2)
	dx := x2 - x1
	dy := y2 - y1
	dz := z2 - z1
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// Elevation angle from radar to target (deg), using ENU
func elevationAngleDeg(radarLat, radarLon, radarH, tgtLat, tgtLon, tgtH float64) float64 {
	xr, yr, zr := llhToEcef(radarLat, radarLon, radarH)
	xt, yt, zt := llhToEcef(tgtLat, tgtLon, tgtH)
	dx := xt - xr
	dy := yt - yr
	dz := zt - zr
	e, n, u := ecefToEnu(dx, dy, dz, radarLat, radarLon)
	horizDist := math.Sqrt(e*e + n*n)
	ang := math.Atan2(u, horizDist)
	return rad2deg(ang)
}

// ------------------ Orizzonte geometrico ------------------

// horizonDistanceSurface returns approximate combined horizon distance (meters)
// for observer height h1 and object height h2 using effective Earth radius factor k
func horizonDistanceSurface(h1, h2, k float64) float64 {
	Re := WGS84_A
	RePrime := k * Re
	return math.Sqrt(2*RePrime*h1) + math.Sqrt(2*RePrime*h2)
}

// ------------------ Propagazione: attenuazione semplice ------------------

// simpleAtmosphericAttenuationDBPerKm returns a crude specific attenuation (dB/km)
// This is a placeholder: for accurate modelling use ITU-R P.676/P.838 tables.
func simpleAtmosphericAttenuationDBPerKm(freqHz, tempC, rh, pressureHPa float64) float64 {
	// crude model: attenuation increases with frequency and humidity
	fGHz := freqHz / 1e9
	// baseline: oxygen + water vapor simplified scaling
	base := 0.01 * math.Pow(fGHz, 2.0) // small at UHF, grows with f^2
	// humidity factor (0..1)
	hf := math.Min(math.Max(rh/100.0, 0.0), 1.0)
	return base * (1.0 + 2.0*hf) // dB/km
}

// ------------------ Equazione radar (monostatica) ------------------

func radarEquationPr(Pt, G, lambda, sigma, R, lossLinear float64) float64 {
	// P_r = Pt * G^2 * lambda^2 * sigma / ((4*pi)^3 * R^4 * L)
	num := Pt * G * G * lambda * lambda * sigma
	den := math.Pow(4*math.Pi, 3) * math.Pow(R, 4) * lossLinear
	return num / den
}

// ------------------ Noise and SNR ------------------
func noisePower(Ts, bandwidth float64) float64 {
	return K_B * Ts * bandwidth
}

// ------------------ Coverage calculation per azimuth ------------------

// computeMaxRangeForAzimuth computes the maximum ground-range (m) at which a
// target of given RCS would be detectable along an azimuth, using simple atm attenuation and horizon.
func computeMaxRangeForAzimuth(r *Radar, env Environment, azimuthDeg float64, dr, rLimit float64) float64 {
	// iterate range steps from dr to rLimit
	maxDetect := 0.0
	for rng := dr; rng <= rLimit; rng += dr {
		// derive point at distance r along azimuth
		lat2, lon2 := destinationPoint(r.Lat, r.Lon, azimuthDeg, rng)
		// sample terrain height (placeholder: sea level) - integrate DEM here
		_ = elevationAngleDeg(r.Lat, r.Lon, r.Height, lat2, lon2, 0.0) // not used here, but could be for clutter?
		// check geometric horizon using k-factor
		dHorizon := horizonDistanceSurface(r.Height, 0.0, r.KFactor)
		if rng > dHorizon {
			// beyond geometric horizon; assume not visible
			break
		}
		// compute 3D slant distance assuming target at ground level
		R3D := distance3D(r.Lat, r.Lon, r.Height, lat2, lon2, 0.0)
		// attenuation per km
		alpha := simpleAtmosphericAttenuationDBPerKm(r.Freq, env.TempC, env.RH, env.PressureHPa)
		AdB := alpha * (R3D / 1000.0)
		Latm := math.Pow(10.0, -AdB/10.0)
		// radar equation
		lambda := C / r.Freq
		Pr := radarEquationPr(r.Pt, r.Gain, lambda, r.Sigma, R3D, r.Loss)
		PrAtt := Pr * Latm
		Pn := noisePower(r.Ts, r.Bandwidth)
		SNR := PrAtt / Pn
		if SNR >= r.SNRThreshold {
			maxDetect = rng
		} else {
			// if not detectable at this r, continue to see if further might be detectable? usually no; can break
			// break to speed up
			break
		}
	}
	return maxDetect
}

// buildCoveragePolygon computes r(theta) for theta in 0..360 and returns vertices lat/lon
func buildCoveragePolygon(r *Radar, env Environment, azStepDeg, dr, rLimit float64) [][]float64 {
	vertices := make([][]float64, 0)
	for theta := 0.0; theta < 360.0; theta += azStepDeg {
		rmax := computeMaxRangeForAzimuth(r, env, theta, dr, rLimit)
		// if rmax == 0, set a very small distance to draw something
		if rmax <= 0 {
			rmax = dr
		}
		lat, lon := destinationPoint(r.Lat, r.Lon, theta, rmax)
		vertices = append(vertices, []float64{lon, lat}) // GeoJSON expects [lon, lat]
	}
	// close polygon
	if len(vertices) > 0 {
		vertices = append(vertices, vertices[0])
	}
	return vertices
}

// export polygon to GeoJSON FeatureCollection
func exportCoverageGeoJSON(r Radar, polygon [][]float64, filename string) error {
	feat := GeoJSONFeature{
		Type: "Feature",
		Properties: map[string]interface{}{
			"radar_id": r.ID,
			"freq_Hz":  r.Freq,
			"pt_W":     r.Pt,
		},
		Geometry: GeoJSONGeometry{
			Type:        "Polygon",
			Coordinates: [][][]float64{},
		},
	}
	// Ugly: Go's typing requires conversion
	coords := make([][]float64, 0)
	coords = append(coords, polygon...)
	feat.Geometry.Coordinates = append(feat.Geometry.Coordinates, coords)
	fc := GeoJSONFeatureCollection{
		Type:     "FeatureCollection",
		Features: []GeoJSONFeature{feat},
	}
	b, err := json.MarshalIndent(fc, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(filename, b, 0644)
}

// ------------------ Esempio main ------------------
func main() {
	// Example radar
	radar := Radar{
		ID:           "R1",
		Lat:          45.0,
		Lon:          12.0,
		Height:       50.0,
		Pt:           1000.0,                    // W
		Gain:         math.Pow(10.0, 30.0/10.0), // 30 dBi
		Freq:         5.6e9,                     // 5.6 GHz
		Sigma:        1.0,                       // m^2
		Loss:         1.0,
		Ts:           500.0,
		Bandwidth:    1e6,
		SNRThreshold: 10.0, // linear ~10 (~10 dB)
		KFactor:      4.0 / 3.0,
	}

	env := Environment{TempC: 15.0, RH: 70.0, PressureHPa: 1013.25}

	azStep := 2.0      // degrees
	dr := 1000.0       // meters step
	rLimit := 200000.0 // 200 km

	polygon := buildCoveragePolygon(&radar, env, azStep, dr, rLimit)
	outfn := "radar_coverage.geojson"
	err := exportCoverageGeoJSON(radar, polygon, outfn)
	if err != nil {
		fmt.Println("Errore export GeoJSON:", err)
	} else {
		fmt.Println("GeoJSON di copertura scritto in:", outfn)
	}
}
